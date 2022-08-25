#ifndef DUNE_NONLINOPT_UNCONSTRAINED_HH
#define DUNE_NONLINOPT_UNCONSTRAINED_HH

#ifdef HAVE_EIGEN
#include <Eigen/Dense>
#endif // HAVE_EIGEN

#include<optional>

#ifdef HAVE_DUNE_COMMON
#include<dune/common/parametertree.hh>
#endif // HAVE_DUNE_COMMON

#include<dune/nonlinopt/conjugationfactor.hh>
#include<dune/nonlinopt/linesearchproposal.hh>
#include<dune/nonlinopt/stoppingcriterion.hh>
#include<dune/nonlinopt/linesearch.hh>
#include<dune/nonlinopt/circularbuffer.hh>

namespace Dune {
  namespace NonlinOpt {

    template<typename,typename> class ConjugationFactorBase;
    template<typename,typename> class LinesearchBase;
    template<typename,typename> class InitialProposalBase;
    template<typename,typename> class LinesearchProposalBase;
    template<typename,typename> class StoppingCriterionBase;
    template<typename,typename> class PreconditionerBase;

    /**
     * @brief A quite general class for unconstrained optimization methods
     *
     * This class provides an interface for several line search based
     * optimization schemes, namely quasi-Newton methods, nonlinear Conjugate
     * Gradient methods, and nonlinear GMRES, as well as hybrid methods that
     * combine two or more of these schemes.
     *
     * @tparam Real  type of parameters and other scalars
     * @tparam Point type of parameter vectors
     */
    template<typename Real, typename Point>
      class UnconstrainedOptimization
      {
        std::size_t verbosity;
        std::size_t max_iter;
        std::size_t reset_iter;
        std::size_t extrap_iter;
        bool        use_scaling;
        bool        quasi_newton;
        bool        extrapolation;
        std::size_t storage_limit;
        bool        gmres_skip_linesearch;
        bool        gmres_reset_on_fail;

        std::unique_ptr<ConjugationFactorBase<Real,Point>>  conj_factor;
        std::unique_ptr<LinesearchBase<Real,Point>>         linesearch;
        std::unique_ptr<InitialProposalBase<Real,Point>>    initial_proposal;
        std::unique_ptr<LinesearchProposalBase<Real,Point>> proposal;
        std::unique_ptr<StoppingCriterionBase<Real,Point>>  criterion;
        std::unique_ptr<PreconditionerBase<Real,Point>>     preconditioner;

        std::size_t iteration = 0;

        Point gradient, prec_grad, direction, old_direction;
        std::optional<Point> old_grad, old_prec_grad, gradient_diff;

        Real value, old_value;                    //! function values
        Real alpha, gmres_alpha;                  //! step sizes
        Real derivative, deriv_alpha, grad_norm2,
             old_grad_norm2, old_grad_grad;       //! scalar products

        CircularBuffer<Point>      point_diff, prec_grad_diff;
        CircularBuffer<Real>       bfgs_factors, scales;
        CircularMatrixBuffer<Real> gmres_matrix;

        public:

        /**
         * @brief Constructor
         *
         * This constructor configures an L-BFGS method with a given window size and
         * optional directional scaling, with the latter enabled by default. Other
         * optimization schemes, like nonlinear CG, BFGS in CG mode, or nonlinear GMRES,
         * can be configured by calling the appropriate methods. The same holds for
         * choosing a custom line search, termination criterion, and so on. Setting the
         * verbosity to zero silences the class, increasing it beyond one prints
         * additional information. The default line search criterion leads to a
         * nonmonotone line search, which is typically more efficient and robust, but
         * may allow a temporary increase in the cost function and be less efficient for
         * certain objective functions. See set_linesearchcriterion if you would like
         * to use a monotone line search instead.
         *
         * @param max_iter_          maximum number of performed iterations
         * @param storage_limit_     window size for L-BFGS, N-GMRES, etc.
         * @param verbosity_         controls amount of information that is printed
         * @param abs_norm_tolerance absolute tolerance for default termination criterion
         * @param rel_norm_tolerance relative tolerance for default termination criterion
         *
         * @see set_bfgs
         * @see set_cg
         * @see set_bfgs_cg
         * @see set_gmres
         * @see set_gmres_bfgs
         * @see set_gmres_bfgs_cg
         */
        UnconstrainedOptimization(std::size_t max_iter_ = 10000, std::size_t storage_limit_ = 10,
            std::size_t verbosity_ = 1, Real abs_norm_tolerance = 1e-6, Real rel_norm_tolerance = 1e-12)
          :
            verbosity(verbosity_), max_iter(max_iter_), storage_limit(std::max(storage_limit_,std::size_t{1})),
            point_diff(storage_limit), prec_grad_diff(storage_limit),
            bfgs_factors(storage_limit), scales(storage_limit), gmres_matrix(storage_limit)
        {
          set_initial_proposal<HagerZhangInitialProposal>();
          set_linesearch<HagerZhangLinesearch>();
          set_linesearchcriterion<NonmonotoneRelaxedWolfeConditions>();
          set_stoppingcriterion<GradientMaxNormCriterion>(abs_norm_tolerance,rel_norm_tolerance);

          set_bfgs();
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts configuration options from a Dune::ParameterTree object,
         * and hands this object to subobjects, thereby guaranteeing that
         * the resulting configuration is self-consistent.
         *
         * @param config ParameterTree object with configuration
         */
        UnconstrainedOptimization(const Dune::ParameterTree& config)
          :
            verbosity            (config.get<std::size_t>("optimization.verbosity",    1)),
            max_iter             (config.get<std::size_t>("optimization.max_iter",     10000)),
            reset_iter           (config.get<std::size_t>("optimization.reset_iter",   max_iter)),
            extrap_iter          (config.get<std::size_t>("optimization.extrap_iter",  1)),
            use_scaling          (config.get<bool>       ("optimization.use_scaling",  true)),
            quasi_newton         (config.get<bool>       ("optimization.quasi_newton", true)),
            extrapolation        (config.get<bool>       ("optimization.extrapolation",false)),
            storage_limit        (config.get<std::size_t>("optimization.storage_limit",10)),
            gmres_skip_linesearch(config.get<bool>       ("gmres.skip_linesearch",     false)),
            gmres_reset_on_fail  (config.get<bool>       ("gmres.reset_on_fail",       false)),
            point_diff(storage_limit), prec_grad_diff(storage_limit),
            bfgs_factors(storage_limit), scales(storage_limit), gmres_matrix(storage_limit)
        {
          // configure nonlinear CG
          const std::string conj_name = config.get<std::string>("optimization.conjugation","none");
          if (conj_name == "hestenes-stiefel" || conj_name == "hs")
            conj_factor = std::make_unique<HestenesStiefelBeta<Real,Point>>(config);
          else if (conj_name == "fletcher-reeves" || conj_name == "fr")
            conj_factor = std::make_unique<FletcherReevesBeta<Real,Point>>(config);
          else if (conj_name == "polak-ribiere-polyak" || conj_name == "prp+")
            conj_factor = std::make_unique<PolakRibierePolyakPlusBeta<Real,Point>>(config);
          else if (conj_name == "dai-yuan" || conj_name == "dy")
            conj_factor = std::make_unique<DaiYuanBeta<Real,Point>>(config);
          else if (conj_name == "hager-zhang" || conj_name == "hz")
            conj_factor = std::make_unique<HagerZhangBeta<Real,Point>>(config);
          else if (conj_name != "none")
            DUNE_THROW(Dune::Exception,"conjugation factor not known, use set_conjugation() method");

          // configure initial line search guess
          const std::string initprop_name = config.get<std::string>("linesearch.initialproposal","hager-zhang");
          if (initprop_name == "hager-zhang")
            initial_proposal = std::make_unique<HagerZhangInitialProposal<Real,Point>>();
          else if (initprop_name == "ceres")
            initial_proposal = std::make_unique<CeresInitialProposal<Real,Point>>();
          else if (initprop_name == "simple")
            initial_proposal = std::make_unique<SimpleInitialProposal<Real,Point>>();
          else
            DUNE_THROW(Dune::Exception,"step width proposal not known, use set_initial_proposal() method");

          // configure subsequent line search guesses
          const std::string prop_name = config.get<std::string>("linesearch.proposal","constant");
          if (prop_name == "absolute")
            proposal = std::make_unique<AbsoluteProposal<Real,Point>>(config);
          else if (prop_name == "relative")
            proposal = std::make_unique<RelativeProposal<Real,Point>>(config);
          else if (prop_name == "constant")
            proposal = std::make_unique<ConstantProposal<Real,Point>>(config);
          else if (prop_name == "adaptive")
            proposal = std::make_unique<AdaptiveProposal<Real,Point>>(config);
          else
            DUNE_THROW(Dune::Exception,"step width proposal not known, use set_proposal() method");

          // configure stopping criterion
          const std::string stop_name = config.get<std::string>("optimization.stoppingcriterion","max_norm");
          if (stop_name == "two_norm")
            criterion = std::make_unique<GradientTwoNormCriterion<Real,Point>>(config);
          else if (stop_name == "max_norm")
            criterion = std::make_unique<GradientMaxNormCriterion<Real,Point>>(config);
          else if (stop_name == "value_decrease")
            criterion = std::make_unique<ValueDecreaseCriterion<Real,Point>>(config);
          else if (stop_name == "combined")
            criterion = std::make_unique<CombinedStoppingCriterion<Real,Point>>(config);
          else
            DUNE_THROW(Dune::Exception,"stopping criterion not known, use set_stoppingcriterion() method");

          // configure line search
          const std::string ls_name = config.get<std::string>("optimization.linesearch","hager-zhang");
          if (ls_name == "hager-zhang")
            linesearch = std::make_unique<HagerZhangLinesearch<Real,Point>>(config);
          else if (ls_name == "cubic")
            linesearch = std::make_unique<CubicLinesearch<Real,Point>>(config);
          else if (ls_name == "quadratic")
            linesearch = std::make_unique<QuadraticLinesearch<Real,Point>>(config);
          else if (ls_name == "derivative-quadratic")
            linesearch = std::make_unique<DerivativeQuadraticLinesearch<Real,Point>>(config);
          else if (ls_name == "secant")
            linesearch = std::make_unique<SecantLinesearch<Real,Point>>(config);
          else if (ls_name == "backtracking")
            linesearch = std::make_unique<BacktrackingLinesearch<Real,Point>>(config);
          else
            DUNE_THROW(Dune::Exception,"line search not known, use set_linesearch() method");
        }
#endif // HAVE_DUNE_COMMON

        /**
         * @brief Define the conjugation parameter for nonlinear CG
         *
         * @tparam Conjugation conjugation parameter class template
         * @tparam Args        types of constructor arguments
         *
         * @param  args        constructor arguments to forward
         */
        template<template<typename,typename> class Conjugation, typename... Args>
          void set_conjugation(Args&&... args)
          {
            conj_factor = std::make_unique<Conjugation<Real,Point>>(std::forward<Args>(args)...);
          }

        /**
         * @brief Define the initial line search proposal function
         *
         * This proposal is used for the first iteration and after resets.
         *
         * @tparam Proposal proposal function class template
         * @tparam Args     types of constructor arguments
         *
         * @param  args     constructor arguments to forward
         */
        template<template<typename,typename> class Proposal, typename... Args>
          void set_initial_proposal(Args&&... args)
          {
            initial_proposal = std::make_unique<Proposal<Real,Point>>(std::forward<Args>(args)...);
          }

        /**
         * @brief Define the line search proposal function
         *
         * This proposal is used for subsequent iterations.
         *
         * @tparam Proposal proposal function class template
         * @tparam Args     types of constructor arguments
         *
         * @param  args     constructor arguments to forward
         */
        template<template<typename,typename> class Proposal, typename... Args>
          void set_proposal(Args&&... args)
          {
            proposal = std::make_unique<Proposal<Real,Point>>(std::forward<Args>(args)...);
          }

        /**
         * @brief Define the line search that should be used
         *
         * @tparam Linesearch line search class template
         * @tparam Args       types of constructor arguments
         *
         * @param  args       constructor arguments to forward
         */
        template<template<typename,typename> class Linesearch, typename... Args>
          void set_linesearch(Args&&... args)
          {
            linesearch = std::make_unique<Linesearch<Real,Point>>(std::forward<Args>(args)...);
          }

        /**
         * @brief Define the line search acceptance criterion
         *
         * @tparam Criterion acceptance criterion class template
         * @tparam Args      types of constructor arguments
         *
         * @param  args      constructor arguments to forward
         */
        template<template<typename,typename> class Criterion, typename... Args>
          void set_linesearchcriterion(Args&&... args)
          {
            linesearch->template set_linesearchcriterion<Criterion>(std::forward<Args>(args)...);
          }

        /**
         * @brief Define the termination criterion
         *
         * @tparam Criterion stopping criterion class template
         * @tparam Args      types of constructor arguments
         *
         * @param  args      constructor arguments to forward
         */
        template<template<typename,typename> class Criterion, typename... Args>
          void set_stoppingcriterion(Args&&... args)
          {
            criterion = std::make_unique<Criterion<Real,Point>>(std::forward<Args>(args)...);
          }

        /**
         * @brief Define a preconditioner that should be applied
         *
         * @tparam Preconditioner preconditioner class template
         * @tparam Args           types of constructor arguments
         *
         * @param  args           constructor arguments to forward
         */
        template<template<typename,typename> class Preconditioner, typename... Args>
          void set_preconditioner(Args&&... args)
          {
            preconditioner = std::make_unique<Preconditioner<Real,Point>>(std::forward<Args>(args)...);
          }

        /**
         * @brief Configure as an L-BFGS method
         *
         * This is the limited-memory Broyden-Fletcher-Goldfarb-Shanno method.
         * It is a generalization of the secant method for multidimensional
         * optimization problems, and creates an approximate inverse Hessian
         * matrix from first-order information about changes in the objective
         * function gradient. By default, we apply approximate eigenvalue
         * scaling based on the Barzilai-Borwein term, which often leads to
         * faster convergence, but may be problematic for certain optimization
         * problems. By default, the method is reset every 1000 iterations.
         *
         * @param use_scaling_ enable / disable directional scaling
         * @param reset_iter_  reset method after this many iterations
         */
        void set_bfgs(bool use_scaling_ = true,
            std::size_t reset_iter_ = std::numeric_limits<std::size_t>::max())
        {
          use_scaling = use_scaling_;
          reset_iter = reset_iter_;

          quasi_newton = true;
          extrapolation = false;
          conj_factor.reset();

          set_proposal<ConstantProposal>();
        }

        /**
         * @brief Configure as a nonlinear CG method
         *
         * This is the nonlinear conjugate gradients method. By default, we
         * use the definition of beta provided by Hager and Zhang. Use
         * set_conjugation to use a different variant of nonlinear CG.
         *
         * @param use_scaling_ enable / disable directional scaling
         * @param reset_iter_  reset method after this many iterations
         *
         * @see set_conjugation
         */
        void set_cg(bool use_scaling_ = false, 
            std::size_t reset_iter_ = std::numeric_limits<std::size_t>::max())
        {
          use_scaling = use_scaling_;
          reset_iter = reset_iter_;

          quasi_newton = false;
          extrapolation = false;
          set_conjugation<HagerZhangBeta>();

          set_proposal<RelativeProposal>(2.);
        }

        /**
         * @brief Configure as an L-BFGS method in CG mode (SCG)
         *
         * This is the L-BFGS method in conjugate gradients mode, i.e., the
         * secant information of the most recent pair of vector differences
         * is not incorporated into the inverse Hessian approximation, and is
         * instead used in a CG-type correction using the previous direction.
         * The resulting method can be interpreted as a nonlinear conjugate
         * gradients method with variable preconditioner matrix (in the form
         * of the BFGS approximate inverse Hessian matrix). By default, the
         * definition of beta by Hestenes and Stiefel is employed.
         *
         * @param use_scaling_ enable / disable directional scaling
         * @param reset_iter_  reset method after this many iterations
         */
        void set_bfgs_cg(bool use_scaling_ = true,
            std::size_t reset_iter_ = std::numeric_limits<std::size_t>::max())
        {
          use_scaling = use_scaling_;
          reset_iter = reset_iter_;

          quasi_newton = true;
          extrapolation = false;
          set_conjugation<HestenesStiefelBeta>();

          set_proposal<RelativeProposal>();
        }

        /**
         * @brief Configure as a nonlinear GMRES method
         *
         * This is the nonlinear generalized minimal residual (N-GMRES)
         * method. It isn't able to generate new search directions by
         * itself, and instead uses another method as "preconditioner" for
         * this purpose, by default a simple steepest descent method. Then,
         * it assembles a local least squares problem based on first-order
         * information about changes in the objective function gradient,
         * and extrapolates into the direction of the minimizer of the
         * local least squares problem. In contrast to the original
         * implementation, we extrapolate from the previous iteration and
         * not the intermediate point if the line search is skipped.
         *
         * @param use_scaling_    enable / disable directional scaling
         * @param reset_iter_     reset method after this many iterations
         * @param extrap_iter_    extrapolate after this many iterations
         * @param reset_on_fail   reset method if direction is no descent direction
         * @param skip_linesearch use fixed step size for preconditioner method
         */
        void set_gmres(bool use_scaling_ = true,
            std::size_t reset_iter_ = std::numeric_limits<std::size_t>::max(),
            std::size_t extrap_iter_ = 1, bool reset_on_fail = false,
            bool skip_linesearch = false)
        {
          use_scaling = use_scaling_;
          reset_iter = reset_iter_;
          extrap_iter = extrap_iter_;
          gmres_skip_linesearch = skip_linesearch;
          gmres_reset_on_fail = reset_on_fail;

          quasi_newton = false;
          extrapolation = true;
          conj_factor.reset();

          set_proposal<RelativeProposal>();
        }

        /**
         * @brief Configure as an N-GMRES method with CG preconditioner
         *
         * This is the nonlinear GMREs method, but with conjugate gradients
         * instead of steepest descent as preconditioner. Using CG as
         * preconditioner often leads to better performance than using CG
         * by itself or N-GMRES with steepest descent. By default, the
         * truncated Polak-Ribière-Polyak (PRP+) variant of nonlinear CG
         * is used, this can be changed using set_conjugation.
         *
         * @param use_scaling_    enable / disable directional scaling
         * @param reset_iter_     reset method after this many iterations
         * @param extrap_iter_    extrapolate after this many iterations
         * @param reset_on_fail   reset method if direction is no descent direction
         * @param skip_linesearch use fixed step size for preconditioner method
         *
         * @see set_gmres
         * @see set_conjugation
         */
        void set_gmres_cg(bool use_scaling_ = true,
            std::size_t reset_iter_ = std::numeric_limits<std::size_t>::max(),
            std::size_t extrap_iter_ = 1, bool reset_on_fail = false,
            bool skip_linesearch = false)
        {
          use_scaling = use_scaling_;
          reset_iter = reset_iter_;
          extrap_iter = extrap_iter_;
          skip_linesearch = skip_linesearch;
          gmres_reset_on_fail = reset_on_fail;

          quasi_newton = false;
          extrapolation = true;
          set_conjugation<PolakRibierePolyakPlusBeta>();

          set_proposal<RelativeProposal>(2.);
        }

        /**
         * @brief Configure as a hybrid L-BFGS / N-GMRES method
         *
         * This configures the solver as a hybrid method that uses the
         * gradient and point differences for both a quasi-Newton
         * L-BFGS step and N-GMRES extrapolation. The frequency of
         * extrapolation steps is reduced compared to stand-alone N-GMRES,
         * since the benefit of extrapolation is less pronounced. Note that
         * this hybrid method is often slower than pure L-BFGS. However,
         * it can be beneficial for a subset of problems.
         *
         * @param use_scaling_    enable / disable directional scaling
         * @param reset_iter_     reset method after this many iterations
         * @param extrap_factor   extrapolate after (storage_limit * extrap_factor) iterations
         * @param reset_on_fail   reset method if direction is no descent direction
         * @param skip_linesearch use fixed step size for preconditioner method
         *
         * @see set_bfgs
         * @see set_gmres
         */
        void set_gmres_bfgs(bool use_scaling_ = true,
            std::size_t reset_iter_ = std::numeric_limits<std::size_t>::max(),
            std::size_t extrap_factor = 1, bool reset_on_fail = false,
            bool skip_linesearch = false)
        {
          use_scaling = use_scaling_;
          reset_iter = reset_iter_;
          extrap_iter = extrap_factor * storage_limit;
          gmres_skip_linesearch = skip_linesearch;
          gmres_reset_on_fail = reset_on_fail;

          quasi_newton = true;
          extrapolation = true;
          conj_factor.reset();

          set_proposal<ConstantProposal>();
        }

        /**
         * @brief Configure as hybrid L-BFGS-CG / N-GMRES method
         *
         * This is the hybrid L-BFGS / N-GMRES method, but with L-BFGS in
         * CG mode (SCG). Again, this combination does not necessarily
         * improve performance compared to a pure L-BFGS method, but may
         * be beneficial for a subset of problems. By default, the
         * truncated Polak-Ribière-Polyak (PRP+) variant of nonlinear CG
         * is used, this can be changed using set_conjugation.
         *
         * @param use_scaling_    enable / disable directional scaling
         * @param reset_iter_     reset method after this many iterations
         * @param extrap_factor   extrapolate after (storage_limit * extrap_factor) iterations
         * @param reset_on_fail   reset method if direction is no descent direction
         * @param skip_linesearch use fixed step size for preconditioner method
         *
         * @see set_gmres_bfgs
         * @see set_conjugation
         */
        void set_gmres_bfgs_cg(bool use_scaling_ = true,
            std::size_t reset_iter_ = std::numeric_limits<std::size_t>::max(),
            std::size_t extrap_factor = 1, bool reset_on_fail = false,
            bool skip_linesearch = false)
        {
          use_scaling = use_scaling_;
          reset_iter = reset_iter_;
          extrap_iter = extrap_factor * storage_limit;
          gmres_skip_linesearch = skip_linesearch;
          gmres_reset_on_fail = reset_on_fail;

          quasi_newton = true;
          extrapolation = true;
          set_conjugation<PolakRibierePolyakPlusBeta>();

          set_proposal<RelativeProposal>();
        }

        /**
         * @brief Provide information about the configuration that is used
         *
         * Prints the names of subalgorithms and termination criteria.
         */
        void report() const
        {
          {
            std::cout << "Method:             ";
            if (extrapolation)
            {
              if (quasi_newton)
              {
                if (conj_factor)
                  std::cout << "hybrid N-GMRES / L-BFGS-CG (SCG)" << std::endl;
                else
                  std::cout << "hybrid N-GMRES / L-BFGS" << std::endl;
              }
              else
              {
                if (conj_factor)
                  std::cout << "nonlinear GMRES (CG preconditioner)" << std::endl;
                else
                  std::cout << "nonlinear GMRES" << std::endl;
              }
            }
            else
            {
              if (quasi_newton)
              {
                if (conj_factor)
                  std::cout << "L-BFGS in CG mode (SCG)" << std::endl;
                else
                  std::cout << "L-BFGS" << std::endl;
              }
              else
              {
                if (conj_factor)
                  std::cout << "nonlinear CG" << std::endl;
                else
                  std::cout << "steepest descent" << std::endl;
              }
            }

            if (conj_factor)
              std::cout << "Conjugation:        " << conj_factor->name() << std::endl;
            else
              std::cout << "Conjugation:        none" << std::endl;
            std::cout << "Initial proposal:   " << initial_proposal->name() << std::endl;
            std::cout << "Line proposal:      " << proposal->name() << std::endl;
            std::cout << "Line search:        " << linesearch->name() << std::endl;
            std::cout << "Line criterion:     " << linesearch->criterion_name() << std::endl;
            std::cout << "Stopping criterion: " << criterion->name() << std::endl;
            if (preconditioner)
              std::cout << "Preconditioner:     " << preconditioner->name() << std::endl;
            else
              std::cout << "Preconditioner:     none" << std::endl;
          }
        }

        /**
         * @brief Access to the current function value
         *
         * @return function value at current point
         */
        Real get_value() const
        {
          return value;
        }

        /**
         * @brief Access to the current gradient or preconditioned gradient
         *
         * @param preconditioned select between gradient and preconditioned gradient
         *
         * @return const reference to gradient or preconditioned gradient
         */
        const Point& get_gradient(bool preconditioned = false) const
        {
          if (preconditioned)
            return prec_grad;

          return gradient;
        }

        /**
         * @brief Reset all information, for starting a new run
         *
         * This function clears all currently stored information,
         * preparing the solver object for a new optimization run
         * with a new starting point, and potentially on a different
         * optimization problem.
         *
         * @param problem problem definition to initialize vectors
         */
        void hard_reset(const ProblemBase<Real,Point>& problem)
        {
          iteration = 0;

          gradient      = problem.zero();
          prec_grad     = problem.zero();
          direction     = problem.zero();
          old_direction = problem.zero();

          if (conj_factor || quasi_newton || extrapolation)
          {
            old_grad = problem.zero();
            old_prec_grad = problem.zero();
          }

          linesearch->hard_reset();
          criterion->hard_reset();

          soft_reset(problem);
        }

        /**
         * @brief Reset information used when determining next direction
         *
         * This function clears all information that is retained from previous
         * iterations, like the vector pairs of BFGS, or the previous direction
         * used in CG. This is used when the method is reset periodically, or
         * when something fails, e.g., an indefinite Hessian approximation in
         * BFGS.
         *
         * @param problem problem definition to initialize vectors
         */
        void soft_reset(const ProblemBase<Real,Point>& problem)
        {
          point_diff.clear();
          prec_grad_diff.clear();
          bfgs_factors.clear();
          gmres_matrix.clear();
          scales.clear();

          old_direction = problem.zero();
        }

        /**
         * @brief Find local minimum of given problem definition
         *
         * @param         problem definition of objective function
         * @param[in,out] point   starting guess, will become local minimizer
         *
         * @return true if termination criterion is fulfilled, else false
         */
        bool apply(const ProblemBase<Real,Point>& problem, Point& point)
        {
          hard_reset(problem);

          bool finished = false;
          while (!finished && iteration < max_iter)
            finished = step(problem,point);

          return finished;
        }

        /**
         * @brief One optimization step
         *
         * @param         problem definition of objective function
         * @param[in,out] point   current iterate, will become next iterate
         *
         * @return true if termination criterion is fulfilled, else false
         */
        bool step(const ProblemBase<Real,Point>& problem, Point& point)
        {
          if (iteration == 0)
          {
            problem.gradient(point,gradient);
            value = problem.value(point,true);

            problem.hook(iteration,point,value,gradient);

            if (criterion->evaluate(old_value,value,gradient))
            {
              if (verbosity >= 1)
                std::cout << "final:     "
                  << " f: " << std::setw(12) << std::scientific << value
                  << " gnorm: " << std::setw(12) << std::scientific << gradient.inf_norm()
                  << std::endl;

              return true;
            }

            prec_grad = gradient;
            if (preconditioner)
              preconditioner->apply(point,prec_grad);
          }
          else if (iteration % reset_iter == 0)
            soft_reset(problem);

          if (verbosity >= 1)
            std::cout << "iter: " << std::setw(5) << iteration
              << " f: " << std::setw(12) << std::scientific << value
              << " gnorm: " << std::setw(12) << std::scientific << gradient.inf_norm()
              << std::endl;

          compute_direction(problem);

          if (iteration % reset_iter == 0)
            alpha = initial_proposal->apply(point,value,gradient);
          else
            alpha = proposal->apply(direction,old_direction,alpha);

          old_value = value;
          old_direction = direction;
          if (conj_factor || quasi_newton || extrapolation || use_scaling)
            old_grad = gradient;

          if (derivative >= 0.)
          {
            if (verbosity >= 2)
              std::cout << "WARNING: not a descent direction, resetting method" << std::endl;

            soft_reset(problem);

            direction = prec_grad;
            direction *= -1.;
            derivative = gradient * direction;

            alpha = initial_proposal->apply(point,value,gradient);
          }

          if (extrapolation && gmres_skip_linesearch)
          {
            const Real dir_norm = std::sqrt(direction * direction);
            alpha = std::min(1e-4/dir_norm,1.);
            Point test_point = point;
            test_point.axpy(direction,alpha);
            Point test_grad;
            problem.gradient(test_point,test_grad);

            point_diff.shift();
            prec_grad_diff.shift();
            bfgs_factors.shift();
            scales.shift();

            point_diff[-1] = point;
            point_diff[-1] -= test_point;

            gradient_diff = gradient;
            *gradient_diff -= test_grad;

            if (preconditioner)
              preconditioner->apply(test_point,test_grad);

            prec_grad_diff[-1] = prec_grad;
            prec_grad_diff[-1] -= test_grad;

            const Real scalar_prod = point_diff[-1] * *gradient_diff;
            bfgs_factors[-1] = 1./scalar_prod;
            if (bfgs_factors[-1] >= 0.)
              scales[-1] = scalar_prod / (prec_grad_diff[-1] * *gradient_diff);
            else
              scales[-1] = 1.;
          }
          else
          {
            try
            {
              linesearch->apply(problem,direction,derivative,point,value,gradient,alpha,deriv_alpha);
            }
            catch(LinesearchFailed const&)
            {
              if (verbosity >= 2)
                std::cout << "WARNING: linesearch failed, resetting method" << std::endl;

              soft_reset(problem);

              direction = prec_grad;
              direction *= -1.;
              derivative = gradient * direction;

              alpha = initial_proposal->apply(point,value,gradient);

              linesearch->apply(problem,direction,derivative,point,value,gradient,alpha,deriv_alpha);
            }
          }

          iteration++;
          problem.hook(iteration,point,value,gradient);

          if (criterion->evaluate(old_value,value,gradient))
          {
            if (verbosity >= 1)
              std::cout << "final:     "
                << " f: " << std::setw(12) << std::scientific << value
                << " gnorm: " << std::setw(12) << std::scientific << gradient.inf_norm()
                << std::endl;

            return true;
          }

          if (conj_factor || quasi_newton || extrapolation || use_scaling)
            old_prec_grad = prec_grad;

          prec_grad = gradient;
          if (preconditioner)
            preconditioner->apply(point,prec_grad);

          if ((!extrapolation && (quasi_newton || use_scaling))
              || (extrapolation && !gmres_skip_linesearch))
          {
            point_diff.shift();
            prec_grad_diff.shift();
            bfgs_factors.shift();
            scales.shift();

            point_diff[-1] = direction;
            point_diff[-1] *= alpha;

            prec_grad_diff[-1] = prec_grad;
            prec_grad_diff[-1] -= *old_prec_grad;

            gradient_diff = gradient;
            *gradient_diff -= *old_grad;

            const Real scalar_prod = point_diff[-1] * *gradient_diff;
            bfgs_factors[-1] = 1./scalar_prod;
            if (bfgs_factors[-1] >= 0.)
              scales[-1] = scalar_prod / (prec_grad_diff[-1] * *gradient_diff);
            else
              scales[-1] = 1.;
          }

          if (extrapolation)
            return extrapolate(problem,point);
          else
            return false;
        }

        private:

        /**
         * @brief Compute search direction
         *
         * This method computes the search direction for the next
         * linesearch. It starts with the negative of the gradient as
         * descent direction, resp. the negative of the preconditioned
         * gradient if a preconditioner is used, and multiplies it with
         * it the L-BFGS approximate inverse Hessian and / or performs
         * a CG directional correction step, depending on configuration.
         *
         * @param problem problem definition, needed for soft reset
         */
        void compute_direction(const ProblemBase<Real,Point>& problem)
        {
          direction = prec_grad;
          direction *= -1.;

          if (!bfgs_factors.empty() && bfgs_factors[-1] < 0.)
          {
            if (verbosity >= 2)
              std::cout << "WARNING: clearing due to negative factor" << std::endl;

            soft_reset(problem);
          }

          if (quasi_newton)
          {
            // skip newest pair if CG mode is used
            // (their information is contained in beta instead)
            const int skip = conj_factor ? 1 : 0;

            std::deque<Real> alphaCoeffs, betaCoeffs;
            for (int i = point_diff.size() - skip; i > 0; i--)
            {
              alphaCoeffs.push_front(bfgs_factors[i-1] * (point_diff[i-1] * direction));
              direction.axpy(prec_grad_diff[i-1],-alphaCoeffs.front());
            }

            if (use_scaling && !scales.empty())
              direction *= scales[-1];

            for (unsigned int i = 0; i + skip < point_diff.size(); i++)
            {
              betaCoeffs.push_back(bfgs_factors[i] * (prec_grad_diff[i] * direction));
              direction.axpy(point_diff[i],alphaCoeffs[i] - betaCoeffs[i]);
            }
          }
          else if (use_scaling && !scales.empty())
              direction *= scales[-1];

          if (conj_factor)
          {
            // here direction is still simply preconditioned gradient
            // (possibly with quasi-Newton variable preconditioner),
            // so use it for relevant scalar products

            old_grad_norm2 = grad_norm2;
            grad_norm2 = - (gradient * direction);
            old_grad_grad = - (*old_grad * direction);

            Real factor = conj_factor->apply(grad_norm2,old_grad_norm2,
                old_grad_grad,derivative,deriv_alpha);
            if (!std::isfinite(factor))
              factor = 0.;

            direction.axpy(old_direction,factor);
          }
          
          derivative = gradient * direction;
        }

        /**
         * @brief Nonlinear GMRES extrapolation step
         *
         * This method performs an N-GMRES step, i.e., it constructs a
         * local least squares surrogate of the objective function based on
         * lagged position / gradient vector pairs, and extrapolates in the
         * direction of the minimizer of this local problem.
         *
         * @param problem problem definition, needed for soft reset
         * @param point   current position, will be updated through extrapolation
         */
        bool extrapolate(const ProblemBase<Real,Point>& problem, Point& point)
        {
          gmres_matrix.shift();
          if (quasi_newton && use_scaling)
          {
            for (unsigned int i = 0; i < gmres_matrix.size() - 1; i++)
            {
              const Real scalar_prod = scales[i] * (*gradient_diff * prec_grad_diff[i]);
              gmres_matrix(i,-1) = gmres_matrix(-1,i) = scalar_prod;
            }
            gmres_matrix(-1,-1) = scales[-1] * (*gradient_diff * prec_grad_diff[-1]);
          }
          else
          {
            for (unsigned int i = 0; i < gmres_matrix.size() - 1; i++)
            {
              const Real scalar_prod = *gradient_diff * prec_grad_diff[i];
              gmres_matrix(i,-1) = gmres_matrix(-1,i) = scalar_prod;
            }
            gmres_matrix(-1,-1) = *gradient_diff * prec_grad_diff[-1];
          }

          if (iteration % extrap_iter != 0)
            return false;

          if (gmres_skip_linesearch || gmres_matrix.size() > 1)
          {
            std::vector<Real> alphas(gmres_matrix.size(),0.);
            std::vector<Real> betas(gmres_matrix.size());

            if (quasi_newton && use_scaling)
            {
              for (unsigned int i = 0; i < gmres_matrix.size(); i++)
                betas[i] = - scales[i] * (prec_grad_diff[i] * gradient);
            }
            else
            {
              for (unsigned int i = 0; i < gmres_matrix.size(); i++)
                betas[i] = - (prec_grad_diff[i] * gradient);
            }

            {
#ifdef HAVE_EIGEN
              Eigen::MatrixXd A(gmres_matrix.size(),gmres_matrix.size());
              Eigen::VectorXd b(gmres_matrix.size());
              for (std::size_t i = 0; i < gmres_matrix.size(); ++i)
              {
                for (std::size_t j = 0; j < gmres_matrix.size(); ++j)
                  A(i,j) = gmres_matrix(i,j);
                b(i) = betas[i];
              }
              Eigen::VectorXd x = A.fullPivHouseholderQr().solve(b);
              for (std::size_t i = 0; i < gmres_matrix.size(); ++i)
                alphas[i] = x(i);
#else
#if HAVE_DUNE_COMMON
              DUNE_THROW(Dune::Exception,"Nonlinear GMRES requires the Eigen3 package for QR.");
#else
              throw std::logic_error("Nonlinear GMRES requires the Eigen3 package for QR.");
#endif // HAVE_DUNE_COMMON
#endif // HAVE_EIGEN
            }

            Point gmres_direction = problem.zero();
            for (std::size_t i = 0; i < gmres_matrix.size(); i++)
              gmres_direction.axpy(point_diff[i],alphas[i]);

            gmres_alpha = 1.;

            derivative = gradient * gmres_direction;

            if (derivative >= 0.)
            {
              if (gmres_skip_linesearch)
              {
                if (verbosity >= 2)
                  std::cout << "WARNING: no descent direction, using original direction" << std::endl;

                gmres_direction = direction;
                derivative = gradient * gmres_direction;

                gmres_alpha = proposal->apply(direction,old_direction,alpha);
              }
              else
              {
                if (gmres_reset_on_fail)
                {
                  if (verbosity >= 2)
                    std::cout << "WARNING: no descent direction, resetting method" << std::endl;
                  soft_reset(problem);
                }
                else if (verbosity >= 2)
                  std::cout << "WARNING: no descent direction, skipping extrapolation" << std::endl;

                return false;
              }
            }

            direction = gmres_direction;

            linesearch->apply(problem,gmres_direction,derivative,point,value,gradient,gmres_alpha,deriv_alpha);

            problem.hook(iteration,point,value,gradient,true);

            if (criterion->evaluate(old_value,value,gradient))
            {
              if (verbosity >= 1)
                std::cout << "final:     "
                  << " f: " << std::setw(12) << std::scientific << value
                  << " gnorm: " << std::setw(12) << std::scientific << gradient.inf_norm()
                  << std::endl;

              return true;
            }

            prec_grad = gradient;
            if (preconditioner)
              preconditioner->apply(point,prec_grad);

            if (gmres_skip_linesearch)
            {
              point_diff[-1] = gmres_direction;
              point_diff[-1] *= gmres_alpha;
            }
            else
              point_diff[-1].axpy(gmres_direction,gmres_alpha);

            prec_grad_diff[-1] = prec_grad;
            prec_grad_diff[-1] -= *old_prec_grad;

            gradient_diff = gradient;
            *gradient_diff -= *old_grad;

            const Real scalar_prod = point_diff[-1] * *gradient_diff;
            bfgs_factors[-1] = 1./scalar_prod;
            if (bfgs_factors[-1] >= 0.)
              scales[-1] = scalar_prod / (prec_grad_diff[-1] * *gradient_diff);
            else
              scales[-1] = 1.;

            if (quasi_newton && use_scaling)
            {
              for (unsigned int i = 0; i < gmres_matrix.size() - 1; i++)
              {
                const Real scalar_prod = scales[i] * (*gradient_diff * prec_grad_diff[i]);
                gmres_matrix(i,-1) = gmres_matrix(-1,i) = scalar_prod;
              }
              gmres_matrix(-1,-1) = scales[-1] * (*gradient_diff * prec_grad_diff[-1]);
            }
            else
            {
              for (unsigned int i = 0; i < gmres_matrix.size() - 1; i++)
              {
                const Real scalar_prod = *gradient_diff * prec_grad_diff[i];
                gmres_matrix(i,-1) = gmres_matrix(-1,i) = scalar_prod;
              }
              gmres_matrix(-1,-1) = *gradient_diff * prec_grad_diff[-1];
            }
          }

          return false;
        }

      };
  }
}

#endif // DUNE_NONLINOPT_UNCONSTRAINED_HH
