#ifndef DUNE_NONLINOPT_LINESEARCHPROPOSAL_HH
#define DUNE_NONLINOPT_LINESEARCHPROPOSAL_HH

#include<cmath>
#include<deque>

namespace Dune {
  namespace NonlinOpt {

    /**
     * @brief Abstract base class for line search initial values in first iteration
     *
     * Derive from this class to provide a custom initial guess for the first
     * iteration.
     */
    template<typename Real, typename Point>
      struct InitialProposalBase
      {
        /**
         * @brief Unique identifier for the initial proposal
         *
         * @return string containing the name
         */
        virtual std::string name() const = 0;

        /**
         * Provide an initial step width based on minimal information
         *
         * @param point   starting position
         * @param descent negative of gradient at starting position
         * @param value   function value at starting position
         *
         * @return initial step width for first line search
         */
        virtual Real apply(const Point& point, const Real value, const Point& descent) const = 0;

        virtual ~InitialProposalBase(){};
      };



    /**
     * @brief Line search proposal for first iteration from CG_DESCENT
     *
     * This is the initial step width of the CG_DESCENT optimization code. It assumes the
     * objective function is quadratic, and returns the ideal step width for this scenario,
     * multiplied with a small constant \f$\psi_0\f$ to keep the initial guess close to the
     * starting position.
     */
    template<typename Real, typename Point>
      class HagerZhangInitialProposal
      : public InitialProposalBase<Real,Point>
      {
        const Real psi_0;

        public:

        /**
         * @brief Constructor
         *
         * @param psi_0_ value for scaling term \f$\psi_0\f$
         */
        HagerZhangInitialProposal(Real psi_0_ = 0.01)
          : psi_0(psi_0_)
        {
          if (psi_0 <= 0.)
            throw std::invalid_argument("psi_0 has to be positive");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the scaling term \f$\psi_0\f$ from the configuration.
         *
         * @param config configuration object
         */
        HagerZhangInitialProposal(const Dune::ParameterTree& config)
          : psi_0(config.get<Real>("linesearch.hz_psi0",0.01))
        {
          if (psi_0 <= 0.)
            DUNE_THROW(Dune::Exception,"psi_0 has to be positive");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "Hager-Zhang initial proposal";
        }

        Real apply(const Point& point, const Real value, const Point& descent) const override
        {
          const Real point_norm = std::sqrt(point*point);
          if (point_norm == 0.)
          {
            const Real abs_func = std::abs(value);
            if (abs_func == 0.)
              return 1.;
            else
              return psi_0 * abs_func / (descent * descent);
          }
          else
            return psi_0 * point.inf_norm() / descent.inf_norm();
        }
      };



    /**
     * @brief Alternative line search proposal for first iteration
     *
     * This is the initial step width of the ceres-solver optimization package. It normalizes
     * the step direction in the maximum norm, so that the largest coordinate change is one.
     */
    template<typename Real, typename Point>
      struct CeresInitialProposal
      : public InitialProposalBase<Real,Point>
      {
        /**
         * @brief Default constructor
         */
        CeresInitialProposal()
        {}

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * @param config configuration object (ignored)
         */
        CeresInitialProposal(const Dune::ParameterTree& config)
        {}
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "ceres initial proposal";
        }

        Real apply(const Point& point, const Real value, const Point& descent) const override
        {
          return std::min(1.,1./descent.inf_norm());
        }
      };



    /**
     * @brief Simple line search proposal for first iteration
     *
     * This initial step width is based on the assumption that the objective function is
     * quadratic, and that its value at the starting point can be decomposed into two
     * contributions of equal magnitude, f(x_0)/2 for the constant part of f and f(x_0)/2
     * for the quadrat term evaluated at x_0:
     * \f[
     * f(x) = \frac{1}{2} f(x_0) + \frac{1}{2} (x - x^\ast)^T (x - x^\ast)
     * \f]
     * Alternatively, this is the ideal step width
     * for a quadratic function with minimum zero, damped by the factor 0.5:
     * \f[
     * f(x) = \frac{1}{2} (x - x^\ast)^T (x - x^\ast)
     * \f]
     */
    template<typename Real, typename Point>
      struct SimpleInitialProposal
      : public InitialProposalBase<Real,Point>
      {
        /**
         * @brief Default constructor
         */
        SimpleInitialProposal()
        {}

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * @param config configuration object (ignored)
         */
        SimpleInitialProposal(const Dune::ParameterTree& config)
        {}
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "simple initial proposal";
        }

        Real apply(const Point& point, const Real value, const Point& descent) const override
        {
          Real abs_func = std::abs(value);
          if (abs_func == 0.)
            abs_func = 1.;

          return abs_func / (descent * descent);
        }
      };



    /**
     * @brief Abstract base class for line search initial values in subsequent iterations
     *
     * Derive from this class to provide a custom initial guess for all line searches
     * except the first one and those after method resets.
     */
    template<typename Real, typename Point>
      struct LinesearchProposalBase
      {
        /**
         * @brief Unique identifier for the initial proposal
         *
         * @return string containing the name
         */
        virtual std::string name() const = 0;

        /**
         * Provide an initial step width based on information from the previous iteration
         *
         * @param direction     current search direction
         * @param old_direction previous search direction
         * @param alpha         step width from previous line search
         *
         * @return initial guess for next line search
         */
        virtual Real apply(const Point& direction,
            const Point& old_direction, const Real alpha) const = 0;

        virtual ~LinesearchProposalBase(){};
      };

    /** @brief Line search proposal that reuses absolute, or ``physical'', step width
     *
     * This proposal assumes that the Euclidean norm of the update vector will be the
     * same as in the previous step, or that the first-order change in the function is
     * the same, as [Nocedal & Wright, Numerical Optimization] put it.
     */
    template<typename Real, typename Point>
      class AbsoluteProposal
      : public LinesearchProposalBase<Real,Point>
      {
        const Real factor;

        public:

        /** @brief Constructor
         *
         * @param factor_ scaling term for proposal
         */
        AbsoluteProposal(Real factor_ = 1.)
          : factor(factor_)
        {
          if (factor <= 0.)
            throw std::invalid_argument("proposal factor has to be positive");
        }

#if HAVE_DUNE_COMMON
        /** @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the multiplicative factor from the configuration.
         *
         * @param config configuration object
         */
        AbsoluteProposal(const Dune::ParameterTree& config)
          : factor(config.get<Real>("linesearch.proposal_factor",1.))
        {
          if (factor <= 0.)
            DUNE_THROW(Dune::Exception,"proposal factor has to be positive");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "absolute proposal";
        }

        Real apply(const Point& direction,
            const Point& old_direction, const Real alpha) const override
        {
          return factor * alpha * std::sqrt(old_direction * old_direction)
            / std::sqrt(direction * direction);
        }

      };



    /**
     * @brief Line search proposal that suggests previous step width, potentially scaled
     *
     * This proposal assumes that the previous step width is a good candidate for the
     * current line search. Optionally, the proposal may be scaled with a small positive
     * constant. A value of one is appropriate for Quasi-Newton methods, while scaling
     * by a factor of two is a popular choice for Conjugate Gradients schemes.
     */
    template<typename Real, typename Point>
      class RelativeProposal
      : public LinesearchProposalBase<Real,Point>
      {
        const Real factor;

        public:

        /** @brief Constructor
         *
         * @param factor_ scaling term for proposal
         */
        RelativeProposal(Real factor_ = 1.)
          : factor(factor_)
        {
          if (factor <= 0.)
            throw std::invalid_argument("proposal factor has to be positive");
        }

#if HAVE_DUNE_COMMON
        /** @brief Constructor
         *
         * Extracts the multiplicative factor from the configuration.
         *
         * @param config configuration object
         */
        RelativeProposal(const Dune::ParameterTree& config)
          : factor(config.get<Real>("linesearch.proposal_factor",1.))
        {
          if (factor <= 0.)
            DUNE_THROW(Dune::Exception,"proposal factor has to be positive");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "relative proposal";
        }

        Real apply(const Point& direction,
            const Point& old_direction, const Real alpha) const override
        {
          return factor * alpha;
        }

      };



    /**
     * @brief Line search proposal that always returns the same constant step width
     *
     * This proposal returns the same step width for each line search. The default
     * value of one for the multiplicative factor will automatically choose the
     * point resulting from uasi-Newton methods, while a value slightly larger than
     * one may help in preventing expensive bracketing steps during line search.
     */
    template<typename Real, typename Point>
      class ConstantProposal
      : public LinesearchProposalBase<Real,Point>
      {
        const Real factor;

        public:

        /** @brief Constructor
         *
         * @param factor_ scaling term for proposal
         */
        ConstantProposal(Real factor_ = 1.)
          : factor(factor_)
        {
          if (factor <= 0.)
            throw std::invalid_argument("proposal factor has to be positive");
        }

#if HAVE_DUNE_COMMON
        /** @brief Constructor
         *
         * Extracts the multiplicative factor from the configuration.
         *
         * @param config configuration object
         */
        ConstantProposal(const Dune::ParameterTree& config)
          : factor(config.get<Real>("linesearch.proposal_factor",1.))
        {
          if (factor <= 0.)
            DUNE_THROW(Dune::Exception,"proposal factor has to be positive");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "constant proposal";
        }

        Real apply(const Point& direction,
            const Point& old_direction, const Real alpha) const override
        {
          return factor;
        }

      };



    /**
     * @brief Line search proposal suggesting step width that would have performed best in the past
     *
     * This proposal keeps track of the initial step widths suggested by the absolute,
     * relative, and constant proposal, and returns the one that would have performed
     * the best across the last few iterations. The window size is configurable.
     */
    template<typename Real, typename Point>
      class AdaptiveProposal
      : public LinesearchProposalBase<Real,Point>
      {
        const std::size_t storage_limit;

        const AbsoluteProposal<Real,Point> absolute_proposal;
        const RelativeProposal<Real,Point> relative_proposal;
        const ConstantProposal<Real,Point> constant_proposal;

        mutable std::deque<Real> absolute, relative, constant, previous;
        mutable std::size_t absolute_best, relative_best, constant_best;

        public:

        /** @brief Constructor
         *
         * @param storage_limit_ number of previous values to store
         * @param factor         scaling term for proposal
         */
        AdaptiveProposal(std::size_t storage_limit_ = 10, Real factor = 1.)
          : storage_limit(storage_limit_),
          absolute_proposal(factor), relative_proposal(factor), constant_proposal(factor)
        {}

#if HAVE_DUNE_COMMON
        /** @brief Constructor
         *
         * Extracts the storage window size and multiplicative factor from the configuration.
         *
         * @param config configuration object
         */
        AdaptiveProposal(const Dune::ParameterTree& config)
          : storage_limit(config.get<std::size_t>("linesearch.storage_limit",10)),
          absolute_proposal(config), relative_proposal(config), constant_proposal(config)
        {}
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "adaptive proposal";
        }

        Real apply(const Point& direction, const Point& old_direction, const Real alpha) const override
        {
          if (absolute.empty())
          {
            absolute.push_back(alpha);
            relative.push_back(alpha);
            constant.push_back(alpha);
          }
          previous.push_back(alpha);

          if (absolute.size() >= storage_limit)
          {
            const Real absolute_diff = std::abs(absolute.front() - previous.front());
            const Real relative_diff = std::abs(relative.front() - previous.front());
            const Real constant_diff = std::abs(constant.front() - previous.front());
            if (constant_diff < relative_diff)
            {
              if (absolute_diff < constant_diff)
                absolute_best--;
              else
                constant_best--;
            }
            else
            {
              if (absolute_diff < relative_diff)
                absolute_best--;
              else
                relative_best--;
            }
            absolute.pop_front();
            relative.pop_front();
            constant.pop_front();
            previous.pop_front();
          }

          const Real absolute_diff = std::abs(absolute.back() - previous.back());
          const Real relative_diff = std::abs(relative.back() - previous.back());
          const Real constant_diff = std::abs(constant.back() - previous.back());
          if (constant_diff < relative_diff)
          {
            if (absolute_diff < constant_diff)
              absolute_best++;
            else
              constant_best++;
          }
          else
          {
            if (absolute_diff < relative_diff)
              absolute_best++;
            else
              relative_best++;
          }

          const Real absolute_prop = absolute_proposal.apply(direction,old_direction,alpha);
          const Real relative_prop = relative_proposal.apply(direction,old_direction,alpha);
          const Real constant_prop = constant_proposal.apply(direction,old_direction,alpha);
          absolute.push_back(absolute_prop);
          relative.push_back(relative_prop);
          constant.push_back(constant_prop);

          if (absolute_best > relative_best && absolute_best > constant_best)
            return absolute_prop;
          else if (constant_best > relative_best)
            return constant_prop;
          else
            return relative_prop;
        }
      };

  }
}

#endif // DUNE_NONLINOPT_LINESEARCHPROPOSAL_HH
