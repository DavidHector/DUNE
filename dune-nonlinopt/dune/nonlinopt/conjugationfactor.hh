#ifndef DUNE_NONLINOPT_CONJUGATIONFACTOR_HH
#define DUNE_NONLINOPT_CONJUGATIONFACTOR_HH

namespace Dune {
  namespace NonlinOpt {

    /**
     * @brief Abstract base class for conjugation factors
     *
     * Derive from this class to provide a custom nonlinear CG method.
     */
    template<typename Real, typename Point>
      struct ConjugationFactorBase
      {
        /**
         * @brief Unique identifier for the conjugation factor
         *
         * @return string containing the name
         */
        virtual std::string name() const = 0;

        /**
         * @brief Compute the conjugation factor
         *
         * @param grad_norm2     square of Euclidean gradient norm
         * @param old_grad_norm2 square of previous Euclidean gradient norm
         * @param old_grad_grad  scalar product of gradient and previous gradient
         * @param derivative     scalar product of gradient with search direction
         * @param deriv_alpha    scalar product of gradient with previous direction
         */
        virtual Real apply(Real grad_norm2, Real old_grad_norm2, Real old_grad_grad,
            Real derivative, Real deriv_alpha) const = 0;

        virtual ~ConjugationFactorBase(){};
      };



    /**
     * @brief Hestenes-Stiefel parameter beta for nonlinear CG
     *
     */
    template<typename Real, typename Point>
      struct HestenesStiefelBeta
      : public ConjugationFactorBase<Real,Point>
      {
        /**
         * @brief Default constructor
         */
        HestenesStiefelBeta()
        {}

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * @param config ParameterTree configuration (ignored)
         */
        HestenesStiefelBeta(const Dune::ParameterTree& config)
        {}
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "Hestenes-Stiefel";
        }

        Real apply(Real grad_norm2, Real old_grad_norm2, Real old_grad_grad,
            Real derivative, Real deriv_alpha) const override
        {
          return (grad_norm2 - old_grad_grad) / (deriv_alpha - derivative);
        }
      };



    /**
     * @brief Fletcher-Reeves parameter beta for nonlinear CG
     *
     */
    template<typename Real, typename Point>
      struct FletcherReevesBeta
      : public ConjugationFactorBase<Real,Point>
      {
        /**
         * @brief Default constructor
         */
        FletcherReevesBeta()
        {}

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * @param config ParameterTree configuration (ignored)
         */
        FletcherReevesBeta(const Dune::ParameterTree& config)
        {}
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "Fletcher-Reeves";
        }

        Real apply(Real grad_norm2, Real old_grad_norm2, Real old_grad_grad,
            Real derivative, Real deriv_alpha) const override
        {
          return grad_norm2 / old_grad_norm2;
        }
      };



    /**
     * @brief Polak-Ribiere-Polyak parameter beta for nonlinear CG
     */
    template<typename Real, typename Point>
      struct PolakRibierePolyakPlusBeta
      : public ConjugationFactorBase<Real,Point>
      {
        /**
         * @brief Default constructor
         */
        PolakRibierePolyakPlusBeta()
        {}

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * @param config ParameterTree configuration (ignored)
         */
        PolakRibierePolyakPlusBeta(const Dune::ParameterTree& config)
        {}
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "Polak-Ribiere-Polyak (PRP+)";
        }

        Real apply(Real grad_norm2, Real old_grad_norm2, Real old_grad_grad,
            Real derivative, Real deriv_alpha) const override
        {
          return std::max(0.,(grad_norm2 - old_grad_grad) / old_grad_norm2);
        }
      };



    /**
     * @brief Dai-Yuan parameter beta for nonlinear CG
     */
    template<typename Real, typename Point>
      struct DaiYuanBeta
      : public ConjugationFactorBase<Real,Point>
      {
        /**
         * @brief Default constructor
         */
        DaiYuanBeta()
        {}

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * @param config ParameterTree configuration (ignored)
         */
        DaiYuanBeta(const Dune::ParameterTree& config)
        {}
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "Dai-Yuan";
        }

        Real apply(Real grad_norm2, Real old_grad_norm2, Real old_grad_grad,
            Real derivative, Real deriv_alpha) const override
        {
          return grad_norm2 / (deriv_alpha - derivative);
        }
      };



    /**
     * @brief Hager-Zhang parameter beta for nonlinear CG
     *
     * This is the nonlinear CG parameter from CG_DESCENT. Note that we
     * choose a different truncation than in the original publication:
     * the original bounding term can't be evaluated if a preconditioner
     * is used. We use the scalar product of gradient and
     * preconditioned gradient as denominator instead.
     */
    template<typename Real, typename Point>
      class HagerZhangBeta
      : public ConjugationFactorBase<Real,Point>
      {
        const Real beta_bound;

        public:

        /**
         * @brief Constructor
         *
         * @param beta_bound_ coefficient for lower cutoff bound
         */
        HagerZhangBeta(Real beta_bound_ = 2.)
          : beta_bound(beta_bound_)
        {
          if (beta_bound <= 0.)
            throw std::invalid_argument("beta bound must be positive");
        }

#if HAVE_DUNE_COMMON
        /**
         * @brief Constructor based on Dune::ParameterTree
         *
         * Extracts the lower bound coefficient from the configuration.
         *
         * @param config ParameterTree configuration (ignored)
         */
        HagerZhangBeta(const Dune::ParameterTree& config)
          : beta_bound(config.get<Real>("optimization.beta_bound",2.))
        {
          if (beta_bound <= 0.)
            DUNE_THROW(Dune::Exception,"beta bound must be positive");
        }
#endif // HAVE_DUNE_COMMON

        std::string name() const override
        {
          return "Hager-Zhang";
        }

        Real apply(Real grad_norm2, Real old_grad_norm2, Real old_grad_grad,
            Real derivative, Real deriv_alpha) const override
        {
          const Real beta_N = ((grad_norm2 - old_grad_grad)
              - 2. * (grad_norm2 - 2.*old_grad_grad + old_grad_norm2)
              * deriv_alpha / (deriv_alpha - derivative)) / (deriv_alpha - derivative);

          const Real bound = beta_bound * derivative / grad_norm2;

          return std::max(beta_N,std::min(0.,bound));
        }
      };

  }
}

#endif // DUNE_NONLINOPT_CONJUGATIONFACTOR_HH
