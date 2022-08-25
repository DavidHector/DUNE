#ifndef DUNE_NONLINOPT_STOPPINGCRITERION_HH
#define DUNE_NONLINOPT_STOPPINGCRITERION_HH

#if HAVE_DUNE_COMMON
#include<dune/common/parametertree.hh>
#endif // HAVE_DUNE_COMMON

namespace Dune {
  namespace NonlinOpt {

    /**
     * @brief Abstract base class for optimization termination criteria
     *
     * Derive from this class to provide a custom termination criterion.
     */
    template<typename Real, typename Point>
      struct StoppingCriterionBase
      {
        /**
         * @brief Unique identifier for the termination criterion
         *
         * @return string containing the name
         */
        virtual std::string name() const = 0;

        /**
         * @brief Reset the internal state of the criterion
         *
         * Criteria that don't have internal state that needs to be
         * reset after a restart may simply implement an empty function.
         */
        virtual void hard_reset() = 0;

        /**
         * @brief Check termination criterion
         *
         * @param old_value previous function value
         * @param value     current function value
         * @param gradient  gradient at current position
         *
         * @return true if criterion is satisfied, else false
         */
        virtual bool evaluate(Real old_value, Real value, const Point& gradient) const = 0;

        virtual ~StoppingCriterionBase(){};
      };



    /**
     * @brief Termination criterion based on Euclidean norm of gradient
     *
     * This criterion stops the optimization method if the Euclidean
     * norm of the gradient falls below a given absolute tolerance,
     * or if the norm falls below a given fraction of the initial
     * gradient norm. The latter test can be disabled by setting the
     * relative tolerance to zero.
     */
    template<typename Real, typename Point>
      class GradientTwoNormCriterion
      : public StoppingCriterionBase<Real, Point>
      {
        const Real abs_norm_tolerance;
        const Real rel_norm_tolerance;

        mutable Real tolerance = 0.;

        public:

        /**
         * @brief Constructor
         *
         * @param abs_norm_tolerance_ absolute tolerance for gradient norm
         * @param rel_norm_tolerance_ tolerance relative to initial norm
         */
        GradientTwoNormCriterion(Real abs_norm_tolerance_ = 1e-6, Real rel_norm_tolerance_ = 1e-12)
          :
            abs_norm_tolerance(abs_norm_tolerance_),
            rel_norm_tolerance(rel_norm_tolerance_)
        {
          if (abs_norm_tolerance <= 0.)
            throw std::invalid_argument("absolute tolerance must be positive");
          if (rel_norm_tolerance < 0.)
            throw std::invalid_argument("relative tolerance must be non-negative");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the norm tolerances from the configuration object.
         *
         * @param config ParameterTree object containing the configuration
         */
        GradientTwoNormCriterion(const Dune::ParameterTree& config)
          :
            abs_norm_tolerance(config.get<Real>("stoppingcriterion.abs_norm_tolerance",1e-6)),
            rel_norm_tolerance(config.get<Real>("stoppingcriterion.rel_norm_tolerance",1e-12))
        {
          if (abs_norm_tolerance <= 0.)
            DUNE_THROW(Dune::Exception,"absolute tolerance must be positive");
          if (rel_norm_tolerance < 0.)
            DUNE_THROW(Dune::Exception,"relative tolerance must be non-negative");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "Euclidean gradient norm";
        }

        void hard_reset() override
        {
          tolerance = 0.;
        }

        bool evaluate(Real old_value, Real value, const Point& gradient) const override
        {
          const Real norm = std::sqrt(gradient * gradient);

          if (tolerance == 0.)
            tolerance = std::max(abs_norm_tolerance, norm * rel_norm_tolerance);

          return (norm <= tolerance);
        }
      };



    /**
     * @brief Termination criterion based on maximum norm of gradient
     *
     * This criterion stops the optimization method if the maximum
     * norm of the gradient falls below a given absolute tolerance,
     * or if the norm falls below a given fraction of the initial
     * gradient norm. The latter test can be disabled by setting the
     * relative tolerance to zero.
     */
    template<typename Real, typename Point>
      class GradientMaxNormCriterion
      : public StoppingCriterionBase<Real,Point>
      {
        const Real abs_norm_tolerance;
        const Real rel_norm_tolerance;

        mutable Real tolerance = 0.;

        public:

        /**
         * @brief Constructor
         *
         * @param abs_norm_tolerance_ absolute tolerance for gradient norm
         * @param rel_norm_tolerance_ tolerance relative to initial norm
         */
        GradientMaxNormCriterion(Real abs_norm_tolerance_ = 1e-6, Real rel_norm_tolerance_ = 1e-12)
          :
            abs_norm_tolerance(abs_norm_tolerance_),
            rel_norm_tolerance(rel_norm_tolerance_)
        {
          if (abs_norm_tolerance <= 0.)
            throw std::invalid_argument("absolute tolerance must be positive");
          if (rel_norm_tolerance < 0.)
            throw std::invalid_argument("relative tolerance must be non-negative");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the norm tolerances from the configuration object.
         *
         * @param config ParameterTree object containing the configuration
         */
        GradientMaxNormCriterion(const Dune::ParameterTree& config)
          :
            abs_norm_tolerance(config.get<Real>("stoppingcriterion.abs_norm_tolerance",1e-6)),
            rel_norm_tolerance(config.get<Real>("stoppingcriterion.rel_norm_tolerance",1e-12))
        {
          if (abs_norm_tolerance <= 0.)
            DUNE_THROW(Dune::Exception,"absolute tolerance must be positive");
          if (rel_norm_tolerance < 0.)
            DUNE_THROW(Dune::Exception,"relative tolerance must be non-negative");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "maximum norm of gradient";
        }

        void hard_reset() override
        {
          tolerance = 0.;
        }

        bool evaluate(Real old_value, Real value, const Point& gradient) const override
        {
          const Real norm = gradient.inf_norm();

          if (tolerance == 0.)
            tolerance = std::max(abs_norm_tolerance, norm * rel_norm_tolerance);

          return (norm <= tolerance);
        }
      };



    /**
     * @brief Termination criterion based on function value decrease
     *
     * This criterion stops the optimization method if there is no
     * further progress, i.e., if the function values don't decrease
     * from iteration to iteration. This criterion may trigger in
     * situations where the method is still far away from a local
     * minimum, but can be useful when stalling optimization runs
     * should be detected.
     */
    template<typename Real, typename Point>
      class ValueDecreaseCriterion
      : public StoppingCriterionBase<Real,Point>
      {
        const Real decrease_tolerance;

        public:

        /**
         * @brief Constructor
         *
         * @param decrease_tolerance_ tolerance for function value decrease
         */
        ValueDecreaseCriterion(Real decrease_tolerance_ = 1e-12)
          : decrease_tolerance(decrease_tolerance_)
        {
          if (decrease_tolerance <= 0.)
            throw std::invalid_argument("decrease tolerance must be positive");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the norm tolerances from the configuration object.
         *
         * @param config ParameterTree object containing the configuration
         */
        ValueDecreaseCriterion(const Dune::ParameterTree& config)
          : decrease_tolerance(config.template get<Real>("stoppingcriterion.decrease_tolerance",1e-12))
        {
          if (decrease_tolerance <= 0.)
            DUNE_THROW(Dune::Exception,"decrease tolerance must be positive");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "function value decrease";
        }

        void hard_reset() override
        {
          // nothing to do, no internal state
        }

        bool evaluate(Real old_value, Real value, const Point& gradient) const override
        {
          return (old_value - value <= (1. + old_value) * decrease_tolerance);
        }
      };



    /**
     * @brief Criterion combining two separate criteria
     *
     * This criterion stops the optimization method if either
     * the Euclidean norm of the gradient is below a given tolerance
     * or the method appears to be stalling. Other combinations of
     * criteria can be created in a similar fashion.
     */
    template<typename Real, typename Point>
      class CombinedStoppingCriterion
      : public StoppingCriterionBase<Real,Point>
      {
        GradientTwoNormCriterion<Real,Point> grad_crit;
        ValueDecreaseCriterion<Real,Point>   val_crit;

        public:

        /**
         * @brief Constructor
         *
         * @param abs_norm_tolerance absolute tolerance for gradient norm
         * @param rel_norm_tolerance tolerance relative to initial norm
         * @param decrease_tolerance tolerance for function value decrease
         */
        CombinedStoppingCriterion(Real abs_norm_tolerance = 1e-6,
            Real rel_norm_tolerance = 1e-12, Real decrease_tolerance = 1e-12)
          : grad_crit(abs_norm_tolerance,rel_norm_tolerance), val_crit(decrease_tolerance)
        {}

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the tolerances from the configuration object.
         *
         * @param config ParameterTree object containing the configuration
         */
        CombinedStoppingCriterion(const Dune::ParameterTree& config)
          : grad_crit(config), val_crit(config)
        {}
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "combined (" + grad_crit.name() + " / " + val_crit.name() + ")";
        }

        void hard_reset() override
        {
          grad_crit.hard_reset();
        }

        bool evaluate(Real old_value, Real value, const Point& gradient) const override
        {
          return (grad_crit.evaluate(old_value,value,gradient)
              || val_crit.evaluate(old_value,value,gradient));
        }
      };

  }
}

#endif // DUNE_NONLINOPT_STOPPINGCRITERION_HH
