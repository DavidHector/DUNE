#ifndef DUNE_NONLINOPT_LINESEARCHCRITERION_HH
#define DUNE_NONLINOPT_LINESEARCHCRITERION_HH

#include<deque>

#if HAVE_DUNE_COMMON
#include<dune/common/parametertree.hh>
#endif // HAVE_DUNE_COMMON

namespace Dune {
  namespace NonlinOpt {

    /**
     * @brief Abstract base class for line search termination criteria
     *
     * Derive from this class to provide a custom termination criterion.
     */
    template<typename Real, typename Point>
      struct LinesearchCriterionBase
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
         * @brief Check whether criterion is fulfilled
         *
         * Determines whether the current proposed step width \f$\alpha\f$
         * produces an acceptable step, typically by checking for
         * things like sufficient function value decrease and curvature
         * conditions. The derivatives are directional derivatives along
         * the search direction, i.e., given as \f$ g^T d\f$, where
         * \f$g\f$ is the gradient in a given point, and \f$d\f$ the
         * current search direction.
         *
         * @param old_value  function value at current point
         * @param derivative derivative along line at current point
         * @param alpha      proposed step size
         * @param test_value function value at proposed point
         * @param test_deriv derivative along line at proposed point
         *
         * @return true if criterion is fulfilled, else false
         */
        virtual bool check(Real old_value, Real derivative, Real alpha,
            Real test_value, Real test_deriv) const
        {
          return check_value(old_value,derivative,alpha,test_value)
            && check_derivative(old_value,derivative,alpha,test_deriv);
        }

        /**
         * @brief Check value part of criterion
         *
         * Ignores the derivative-based part of the criterion (if any).
         * This allows backing out early if the value-based part
         * rejects the point, without calculating the gradient. The
         * costs for gradient computations are often significant, and
         * this avoids computations that will be discarded anyway.
         *
         * @param old_value  function value at current point
         * @param derivative derivative along line at current point
         * @param alpha      proposed step size
         * @param test_value function value at proposed point
         *
         * @see check
         */
        virtual bool check_value(Real old_value, Real derivative,
            Real alpha, Real test_value) const = 0;

        /**
         * @brief Check derivative part of criterion
         *
         * Ignores the value-based part of the criterion (if any).
         * This allows backing out early if the derivative-based part
         * rejects the point, without calculating the function value.
         * Many problems calculate the function value as a
         * byproduct when calculating the gradient, so this is not
         * that important.
         *
         * @param old_value  function value at current point
         * @param derivative derivative along line at current point
         * @param alpha      proposed step size
         * @param test_deriv derivative along line at proposed point
         *
         * @see check
         */
        virtual bool check_derivative(Real old_value, Real derivative,
            Real alpha, Real test_deriv) const = 0;

        virtual ~LinesearchCriterionBase(){};
      };



    /**
     * @brief Sufficient descent condition (Armijo rule)
     *
     * This condition ensures that the function value decreases by a
     * sufficient amount when jumping to the proposed point, relative
     * to what a linear extrapolation would produce. Note that the
     * Armijo rule can't guarantee that the next direction of a
     * quasi-Newton method will be a descent direction as needed,
     * e.g., by the L-BFGS method, and is therefore a poor choice in
     * that case.
     */
    template<typename Real, typename Point>
      class ArmijoRule
      : public LinesearchCriterionBase<Real,Point>
      {
        const Real armijo_factor;

        public:


        /**
         * @brief Constructor
         *
         * @param armijo_factor_ parameter for sufficient descent condition
         */
        ArmijoRule(Real armijo_factor_ = 0.1)
          : armijo_factor(armijo_factor_)
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            throw std::invalid_argument("Armijo factor has to be in (0,1)");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the Armijo constant from the configuration object.
         *
         * @param config configuration object
         */
        ArmijoRule(const Dune::ParameterTree& config)
          : armijo_factor(config.get<Real>("linesearch.armijo_factor",0.1))
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            DUNE_THROW(Dune::Exception,"Armijo factor has to be in (0,1)");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "Armijo";
        }

        void hard_reset() override
        {
          // no internal state, nothing to do
        }

        bool check_value(Real old_value, Real derivative,
            Real alpha, Real test_value) const override
        {
          if (!std::isfinite(test_value))
            return false;

          const Real armijo_value = old_value + alpha * armijo_factor * derivative;

          return (test_value <= armijo_value);
        }

        bool check_derivative(Real old_value, Real derivative,
            Real alpha, Real test_deriv) const override
        {
          return true;
        }
      };



    /**
     * @brief A nonmonotone variant of the Armijo rule
     *
     * This condition accepts the proposed point if its function
     * value is smaller than the weighted average of the previous
     * iterations. As a consequence, the function value may
     * temporarily increase, and the method becomes nonmonotone.
     */
    template<typename Real, typename Point>
      class NonmonotoneArmijoRule
      : public LinesearchCriterionBase<Real,Point>
      {
        const Real armijo_factor;
        const Real eta;

        mutable Real Q = 1;
        mutable Real C = std::numeric_limits<Real>::max();

        public:

        /**
         * @brief Constructor
         *
         * @param armijo_factor_ parameter for sufficient descent condition
         * @param eta_           weighting factor for average of previous values
         */
        NonmonotoneArmijoRule(Real armijo_factor_ = 0.1, Real eta_ = 1.)
          : armijo_factor(armijo_factor_), eta(eta_)
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            throw std::invalid_argument("Armijo factor has to be in (0,1)");
          if (eta < 0. || eta > 1.)
            throw std::invalid_argument("eta has to be in [0,1]");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the Armijo constant and window size from the
         * configuration object.
         *
         * @param config configuration object
         */
        NonmonotoneArmijoRule(const Dune::ParameterTree& config)
          :
            armijo_factor(config.get<Real>("linesearch.armijo_factor",0.1)),
            eta          (config.get<Real>("linesearch.hz_eta",       1.))
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            DUNE_THROW(Dune::Exception,"Armijo factor has to be in (0,1)");
          if (eta < 0. || eta > 1.)
            DUNE_THROW(Dune::Exception,"eta has to be in [0,1]");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "nonmonotone Armijo";
        }

        void hard_reset() override
        {
          Q = 1;
          C = std::numeric_limits<Real>::max();
        }

        bool check_value(Real old_value, Real derivative,
            Real alpha, Real test_value) const override
        {
          if (!std::isfinite(test_value))
            return false;

          if (Q == 1)
            C = old_value;

          const Real armijo_value = C + alpha * armijo_factor * derivative;

          if (test_value <= armijo_value)
          {
            C = eta * Q * C + test_value;
            Q = eta * Q + 1;
            C /= Q;

            return true;
          }
          else
            return false;
        }

        bool check_derivative(Real old_value, Real derivative,
            Real alpha, Real test_deriv) const override
        {
          return true;
        }
      };



    /**
     * @brief Sufficient descent and curvature conditions (Wolfe conditions)
     *
     * The Wolfe conditions combine the sufficient decrease condition
     * (Armijo rule) with a curvature condition, which guarantees
     * that the directions produced by quasi-Newton methods, e.g.,
     * L-BFGS, are always descent directions. This is the standard
     * criterion for BFGS and similar methods.
     */
    template<typename Real, typename Point>
      class WolfeConditions
      : public LinesearchCriterionBase<Real,Point>
      {
        const Real armijo_factor;
        const Real wolfe_factor;

        public:

        /**
         * @brief Constructor
         *
         * @param armijo_factor_ parameter for sufficient descent condition
         * @param wolfe_factor_  parameter for curvature condition
         */
        WolfeConditions(Real armijo_factor_ = 0.1, Real wolfe_factor_ = 0.9)
          : armijo_factor(armijo_factor_), wolfe_factor(wolfe_factor_)
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            throw std::invalid_argument("Armijo factor has to be in (0,1)");
          if (wolfe_factor <= armijo_factor || wolfe_factor >= 1.)
            throw std::invalid_argument("Wolfe factor has to be in (0,1) and larger than Armijo factor");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the Armijo and Wolfe constants from the configuration object.
         *
         * @param config configuration object
         */
        WolfeConditions(const Dune::ParameterTree& config)
          :
            armijo_factor(config.get<Real>("linesearch.armijo_factor",0.1)),
            wolfe_factor (config.get<Real>("linesearch.wolfe_factor", 0.9))
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            DUNE_THROW(Dune::Exception,"Armijo factor has to be in (0,1)");
          if (wolfe_factor <= armijo_factor || wolfe_factor >= 1.)
            DUNE_THROW(Dune::Exception,"Wolfe factor has to be in (0,1) and larger than Armijo factor");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "Wolfe";
        }

        void hard_reset() override
        {
          // no internal state, nothing to do
        }

        bool check_value(Real old_value, Real derivative,
            Real alpha, Real test_value) const override
        {
          if (!std::isfinite(test_value))
            return false;

          const Real armijo_value = old_value + alpha * armijo_factor * derivative;

          return test_value <= armijo_value;
        }

        bool check_derivative(Real old_value, Real derivative,
            Real alpha, Real test_deriv) const override
        {
          const Real wolfe_value  = wolfe_factor * derivative;

          return test_deriv >= wolfe_value;
        }

      };



    /**
     * @brief Approximation of the Wolfe conditions
     *
     * The sufficient decrease condition may fail due to numerical
     * errors and finite precision, e.g., when near the optimum.
     * Hager and Zhang propose an approximation of the Wolfe
     * conditions that replaces the Armijo rule with a complementary
     * curvature condition, giving a range of allowed values for
     * the directional derivative. This version also ensures that
     * subsequent directions are descent directions.
     *
     * @see WolfeConditions
     */
    template<typename Real, typename Point>
      class ApproximateWolfeConditions
      : public LinesearchCriterionBase<Real,Point>
      {
        const Real armijo_factor;
        const Real wolfe_factor;
        const Real epsilon;

        public:

        /**
         * @brief Constructor
         *
         * @param armijo_factor_ parameter for sufficient descent condition
         * @param wolfe_factor_  parameter for curvature condition
         * @param epsilon_       estimate of function value error
         */
        ApproximateWolfeConditions(Real armijo_factor_ = 0.1,
            Real wolfe_factor_ = 0.9, Real epsilon_ = 1e-6)
          :
            armijo_factor(armijo_factor), wolfe_factor(wolfe_factor_),
            epsilon(epsilon_)
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            throw std::invalid_argument("Armijo factor has to be in (0,1)");
          if (wolfe_factor <= armijo_factor || wolfe_factor >= 1.)
            throw std::invalid_argument("Wolfe factor has to be in (0,1) and larger than Armijo factor");
          if (epsilon < 0.)
            throw std::invalid_argument("estimated function error epsilon has to be non-negative");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the parameters from the configuration object.
         *
         * @param config configuration object
         */
        ApproximateWolfeConditions(const Dune::ParameterTree& config)
          :
            armijo_factor(config.get<Real>("linesearch.armijo_factor",0.1)),
            wolfe_factor (config.get<Real>("linesearch.wolfe_factor", 0.9)),
            epsilon      (config.get<Real>("linesearch.hz_epsilon",   1e-6))
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            DUNE_THROW(Dune::Exception,"Armijo factor has to be in (0,1)");
          if (wolfe_factor <= armijo_factor || wolfe_factor >= 1.)
            DUNE_THROW(Dune::Exception,"Wolfe factor has to be in (0,1) and larger than Armijo factor");
          if (epsilon < 0.)
            DUNE_THROW(Dune::Exception,"estimated function error epsilon has to be non-negative");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "approximate Wolfe";
        }

        void hard_reset() override
        {
          // no internal state, nothing to do
        }

        bool check_value(Real old_value, Real derivative,
            Real alpha, Real test_value) const override
        {
          return test_value <= old_value + epsilon;
        }

        bool check_derivative(Real old_value, Real derivative,
            Real alpha, Real test_deriv) const override
        {
          const Real wolfe_value     = wolfe_factor * derivative;
          const Real alt_wolfe_value = (2.*armijo_factor - 1.) * derivative;

          return test_deriv >= wolfe_value && test_deriv <= alt_wolfe_value;
        }

      };



    /**
     * @brief Combination of Wolfe conditions and their approximation
     *
     * This criterion accepts the proposed point when either the
     * original Wolfe conditions or the approximate Wolfe conditions
     * would accept it.
     *
     * @see WolfeConditions
     * @see ApproximateWolfeConditions
     */
    template<typename Real, typename Point>
      class RelaxedWolfeConditions
      : public LinesearchCriterionBase<Real,Point>
      {
        const Real armijo_factor;
        const Real wolfe_factor;
        const Real epsilon;

        public:

        /**
         * @brief Constructor
         *
         * @param armijo_factor_ parameter for sufficient descent condition
         * @param wolfe_factor_  parameter for curvature condition
         * @param epsilon_       estimate of function value error
         */
        RelaxedWolfeConditions(Real armijo_factor_ = 0.1,
            Real wolfe_factor_ = 0.9, Real epsilon_ = 1e-6)
          :
            armijo_factor(armijo_factor_), wolfe_factor(wolfe_factor_),
            epsilon(epsilon_)
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            throw std::invalid_argument("Armijo factor has to be in (0,1)");
          if (wolfe_factor <= armijo_factor || wolfe_factor >= 1.)
            throw std::invalid_argument("Wolfe factor has to be in (0,1) and larger than Armijo factor");
          if (epsilon < 0.)
            throw std::invalid_argument("estimated function error epsilon has to be non-negative");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the parameters from the configuration object.
         *
         * @param config configuration object
         */
        RelaxedWolfeConditions(const Dune::ParameterTree& config)
          :
            armijo_factor(config.get<Real>("linesearch.armijo_factor",0.1)),
            wolfe_factor (config.get<Real>("linesearch.wolfe_factor", 0.9)),
            epsilon      (config.get<Real>("linesearch.hz_epsilon",   1e-6))
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            DUNE_THROW(Dune::Exception,"Armijo factor has to be in (0,1)");
          if (wolfe_factor <= armijo_factor || wolfe_factor >= 1.)
            DUNE_THROW(Dune::Exception,"Wolfe factor has to be in (0,1) and larger than Armijo factor");
          if (epsilon < 0.)
            DUNE_THROW(Dune::Exception,"estimated function error epsilon has to be non-negative");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "relaxed Wolfe";
        }

        void hard_reset() override
        {
          // no internal state, nothing to do
        }

        bool check(Real old_value, Real derivative,
            Real alpha, Real test_value, Real test_deriv) const override
        {
          if (!std::isfinite(test_value))
            return false;

          const Real armijo_value    = old_value + alpha * armijo_factor * derivative;
          const Real wolfe_value     = wolfe_factor * derivative;
          const Real alt_wolfe_value = (2.*armijo_factor - 1.) * derivative;

          return ((test_value <= old_value + epsilon && test_deriv <= alt_wolfe_value)
              || test_value <= armijo_value) && test_deriv >= wolfe_value;
        }

        bool check_value(Real old_value, Real derivative,
            Real alpha, Real test_value) const override
        {
#if HAVE_DUNE_COMMON
          DUNE_THROW(Dune::Exception,"relaxed Wolfe conditions require derivatives");
#else
          throw std::logic_error("relaxed Wolfe conditions require derivatives");
#endif // HAVE_DUNE_COMMON
        }

        bool check_derivative(Real old_value, Real derivative,
            Real alpha, Real test_deriv) const override
        {
          const Real wolfe_value  = wolfe_factor * derivative;

          return test_deriv >= wolfe_value;
        }

      };



    /**
     * @brief Nonmonotone variant of the relaxed Wolfe conditions
     *
     * This criterion accepts the proposed point if the relaxed
     * Wolfe conditons are fulfilled, with the current function
     * value replaced by a weighted average of the previous function
     * values. As a consequence, the accepted new function value
     * may be larger than the current one, and the resulting method
     * is nonmonotone.
     *
     * @see RelaxedWolfeConditions
     * @see NonmonotoneArmijoRule
     */
    template<typename Real, typename Point>
      class NonmonotoneRelaxedWolfeConditions
      : public LinesearchCriterionBase<Real,Point>
      {
        const Real armijo_factor;
        const Real wolfe_factor;
        const Real epsilon;
        const Real eta;

        mutable Real Q = 1;
        mutable Real C = std::numeric_limits<Real>::max();

        public:

        /**
         * @brief Constructor
         *
         * @param armijo_factor_ parameter for sufficient descent condition
         * @param wolfe_factor_  parameter for curvature condition
         * @param epsilon_       estimate of function value error
         * @param eta_           weighting factor for average of previous values
         */
        NonmonotoneRelaxedWolfeConditions(Real armijo_factor_ = 0.1,
            Real wolfe_factor_ = 0.9, Real epsilon_ = 1e-6, Real eta_ = 1.)
          :
            armijo_factor(armijo_factor_), wolfe_factor(wolfe_factor_),
            epsilon(epsilon_), eta(eta_)
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            throw std::invalid_argument("Armijo factor has to be in (0,1)");
          if (wolfe_factor <= armijo_factor || wolfe_factor >= 1.)
            throw std::invalid_argument("Wolfe factor has to be in (0,1) and larger than Armijo factor");
          if (epsilon < 0.)
            throw std::invalid_argument("estimated function error epsilon has to be non-negative");
          if (eta < 0. || eta > 1.)
            throw std::invalid_argument("eta has to be in [0,1]");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the parameters from the configuration object.
         *
         * @param config configuration object
         */
        NonmonotoneRelaxedWolfeConditions(const Dune::ParameterTree& config)
          :
            armijo_factor(config.get<Real>("linesearch.armijo_factor",0.1)),
            wolfe_factor (config.get<Real>("linesearch.wolfe_factor", 0.9)),
            epsilon      (config.get<Real>("linesearch.hz_epsilon",   1e-6)),
            eta          (config.get<Real>("linesearch.hz_eta",       1.))
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            DUNE_THROW(Dune::Exception,"Armijo factor has to be in (0,1)");
          if (wolfe_factor <= armijo_factor || wolfe_factor >= 1.)
            DUNE_THROW(Dune::Exception,"Wolfe factor has to be in (0,1) and larger than Armijo factor");
          if (epsilon < 0.)
            DUNE_THROW(Dune::Exception,"estimated function error epsilon has to be non-negative");
          if (eta < 0. || eta > 1.)
            DUNE_THROW(Dune::Exception,"eta has to be in [0,1]");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "nonmonotone relaxed Wolfe";
        }

        void hard_reset() override
        {
          Q = 1;
          C = std::numeric_limits<Real>::max();
        }

        bool check(Real old_value, Real derivative,
            Real alpha, Real test_value, Real test_deriv) const override
        {
          if (!std::isfinite(test_value))
            return false;

          if (Q == 1)
            C = old_value;

          const Real armijo_value    = C + alpha * armijo_factor * derivative;
          const Real wolfe_value     = wolfe_factor * derivative;
          const Real alt_wolfe_value = (2.*armijo_factor - 1.) * derivative;

          if (((test_value <= C + epsilon && test_deriv <= alt_wolfe_value)
              || test_value <= armijo_value) && test_deriv >= wolfe_value)
          {
            C = eta * Q * C + test_value;
            Q = eta * Q + 1;
            C /= Q;

            return true;
          }
          else
            return false;
        }

        bool check_value(Real old_value, Real derivative,
            Real alpha, Real test_value) const override
        {
#if HAVE_DUNE_COMMON
          DUNE_THROW(Dune::Exception,"relaxed Wolfe conditions require derivatives");
#else
          throw std::logic_error("relaxed Wolfe conditions require derivatives");
#endif // HAVE_DUNE_COMMON
        }

        bool check_derivative(Real old_value, Real derivative,
            Real alpha, Real test_deriv) const override
        {
          const Real wolfe_value  = wolfe_factor * derivative;

          return test_deriv >= wolfe_value;
        }

      };



    /**
     * @brief Strong Wolfe conditions
     *
     * The strong Wolfe conditions combine the curvature condition
     * of the Wolfe conditions with an additional condition that
     * ensures that the absolute value of the directional derivative
     * at the proposed point is bounded, thereby excluding points
     * that are very far away from stationary points along the line.
     */
    template<typename Real, typename Point>
      class StrongWolfeConditions
      : public LinesearchCriterionBase<Real,Point>
      {
        const Real armijo_factor;
        const Real wolfe_factor;

        public:

        /**
         * @brief Constructor
         *
         * @param armijo_factor_ parameter for sufficient descent condition
         * @param wolfe_factor_  parameter for curvature condition
         */
        StrongWolfeConditions(Real armijo_factor_ = 0.1, Real wolfe_factor_ = 0.9)
          : armijo_factor(armijo_factor_), wolfe_factor(wolfe_factor_)
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            throw std::invalid_argument("Armijo factor has to be in (0,1)");
          if (wolfe_factor <= armijo_factor || wolfe_factor >= 1.)
            throw std::invalid_argument("Wolfe factor has to be in (0,1) and larger than Armijo factor");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the Armijo and Wolfe constants from the configuration object.
         *
         * @param config configuration object
         */
        StrongWolfeConditions(const Dune::ParameterTree& config)
          :
            armijo_factor(config.get<Real>("linesearch.armijo_factor",0.1)),
            wolfe_factor (config.get<Real>("linesearch.wolfe_factor", 0.9))
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            DUNE_THROW(Dune::Exception,"Armijo factor has to be in (0,1)");
          if (wolfe_factor <= armijo_factor || wolfe_factor >= 1.)
            DUNE_THROW(Dune::Exception,"Wolfe factor has to be in (0,1) and larger than Armijo factor");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "strong Wolfe";
        }

        void hard_reset() override
        {
          // no internal state, nothing to do
        }

        bool check_value(Real old_value, Real derivative,
            Real alpha, Real test_value) const override
        {
          if (!std::isfinite(test_value))
            return false;

          const Real armijo_value = old_value + alpha * armijo_factor * derivative;

          return test_value <= armijo_value;
        }

        bool check_derivative(Real old_value, Real derivative,
            Real alpha, Real test_deriv) const override
        {
          const Real wolfe_value  = wolfe_factor * derivative;

          return std::abs(test_deriv) <= std::abs(wolfe_value);
        }
      };



    /**
     * @brief Goldstein Conditions
     *
     * The Goldstein conditions are an alternative to the
     * Wolfe conditions, replacing the curvature condition
     * with a complementary Armijo-type condition that controls
     * the step size from below. Just like the Armijo rule, they
     * can't guarantee that subsequent directions will be descent
     * directions, and therefore the Wolfe conditions are usually
     * preferred for quasi-Newton methods.
     */
    template<typename Real, typename Point>
      class GoldsteinConditions
      : public LinesearchCriterionBase<Real,Point>
      {
        const Real armijo_factor;

        public:

        /**
         * @brief Constructor
         *
         * @param armijo_factor_ parameter for sufficient descent condition
         */
        GoldsteinConditions(Real armijo_factor_ = 0.1)
          : armijo_factor(armijo_factor_)
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            throw std::invalid_argument("Armijo factor has to be in (0,1)");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the Armijo constant from the configuration object.
         *
         * @param config configuration object
         */
        GoldsteinConditions(const Dune::ParameterTree& config)
          : armijo_factor(config.get<Real>("linesearch.armijo_factor",0.1))
        {
          if (armijo_factor <= 0. || armijo_factor >= 1.)
            DUNE_THROW(Dune::Exception,"Armijo factor has to be in (0,1)");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "Goldstein";
        }

        void hard_reset() override
        {
          // no internal state, nothing to do
        }

        bool check_value(Real old_value, Real derivative,
            Real alpha, Real test_value) const override
        {
          if (!std::isfinite(test_value))
            return false;

          const Real armijo_value     = old_value + alpha * armijo_factor * derivative;
          const Real alt_armijo_value = old_value + alpha * (1. - armijo_factor) * derivative;

          return (test_value <= armijo_value && alt_armijo_value <= test_value);
        }

        bool check_derivative(Real old_value, Real derivative,
            Real alpha, Real test_deriv) const override
        {
          return true;
        }
      };

  }
}

#endif // DUNE_NONLINOPT_LINESEARCHCRITERION_HH
