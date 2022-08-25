#ifndef DUNE_NONLINOPT_USERINTERFACE_HH
#define DUNE_NONLINOPT_USERINTERFACE_HH

#ifdef HAVE_MUPARSER
#include "muParser.h"
#endif // HAVE_MUPARSER

#include<dune/nonlinopt/vectorclass.hh>
namespace Dune {
  namespace NonlinOpt {

    /**
     * @brief Abstract base class for optimization problems
     *
     * Derive from this class to define an optimization problem. The
     * optimization code assumes that minimization is the task, i.e.,
     * maximization problems should be reformulated as minimization
     * problems by providing the negative of function values and gradients.
     *
     * @tparam Real  data type of scalars and function values
     * @tparam Point data type of parameter vectors
     */
    template<typename Real, typename Point>
      struct ProblemBase
      {
        /**
         * @brief Evaluate the objective function value
         *
         * This represents the evaluation of the function that should be
         * minimized. Complicated models may have internal state, e.g.,
         * those based on the evaluation of PDEs may contain internal
         * finite element solvers and solutions. The Boolean flag signals
         * that the function value is requested at the same position as
         * the most recent gradient or dir_deriv evaluation, and costly
         * recomputations may be avoided if the function value is still
         * available in the internal state.
         *
         * @param point      vector of function arguments
         * @param subsequent whether arguments are same as previous gradient / line_deriv call
         *
         * @return function value at given position
         *
         * @see dir_deriv
         * @see gradient
         * @see value_and_deriv
         */
        virtual Real value(const Point& point, bool subsequent = false) const = 0;

        /**
         * @brief Evaluate the objective function directional derivative
         *
         * This  represents the evaluation of the directional derivative
         * of the objective function in the given point and given direction.
         * For certain models this produces the function value and gradient as
         * a natural byproduct, e.g., when using adjoint PDEs to compute the
         * sensitivities with regard to a PDE-based model, the model itself
         * has to be evaluated first, and the directional derivative is
         * computed using the gradient. The Boolean flag signals that the
         * derivative is requested at the same position as the most recent
         * function evaluation, and intermediate values that are still
         * available in the internal state may be used during derivative
         * evaluation. The method should update its gradient argument if it
         * is available.
         *
         * @param      point      vector of function arguments
         * @param      direction  direction of directional derivative
         * @param[out] gradient   gradient at given position, if available
         * @param      subsequent whether arguments are same as previous value call
         *
         * @return directional derivative at given position in given direction
         *
         * @see value
         * @see gradient
         * @see value_and_deriv
         */
        virtual Real dir_deriv(const Point& point, const Point& direction,
            Point& gradient, bool subsequent = false) const = 0;

        /**
         * @brief Evaluate the objective function gradient
         *
         * This represents the evaluation of the gradient of the function
         * that should be minimized. The Boolean flag signals that the
         * gradient is requested at the same position as the most recent
         * directional derivative evaluation. For adjoint-based derivation
         * this means that the gradient has already been computed and may be
         * skipped here. For finite differences this means that the directional
         * derivative may be reused when computing the full gradient.
         *
         * @param      point      vector of function arguments
         * @param[out] gradient   gradient at given position
         * @param      direction  direction of directional derivative, if available
         * @param      derivative value of directional derivative, if available
         *
         * @see value
         * @see dir_deriv
         * @see value_and_deriv
         */
        virtual void gradient(const Point& point, Point& gradient,
            const Point* const direction = nullptr, Real derivative = 0.) const = 0;

        /**
         * @brief Evaluate both the objective function and its directional derivative
         *
         * This method evaluates both the function and its directional derivative
         * in one function call, and additionally updates the gradient if it
         * becomes available during computation. Superficially, this is identical
         * to calling both the value and the dir_deriv method. However, there are
         * subtle differences in the case of more advanced line searches, since
         * the line search may back out early if a candidate step width can be
         * rejected based on, say, the function value alone, but the gradient
         * would still be used for decisions if it were available. By default,
         * evaluations are separate, but combined evaluations can be enabled
         * if the line search supports it.
         *
         * @param      point     vector of function arguments
         * @param      direction direction for directional derivative
         * @param[out] gradient  gradient at given position, if available
         *
         * @return pair of function value and directional derivative
         *
         * @see value
         * @see dir_deriv
         * @see gradient
         */
        virtual std::pair<Real,Real> value_and_deriv(const Point& point,
            const Point& direction, Point& gradient) const = 0;

        /**
         * @brief Produce a zero vector of the correct dimension
         *
         * This method returns a vector that represents the origin of the
         * vector space. It is used to initialize vectors during the first
         * iteration to ensure that they have the right dimension.
         *
         * @return vector of parameters with length dim
         *
         * @see dim
         */
        virtual Point zero() const = 0;

        /**
         * @brief Dimension of the optimization problem
         *
         * This is the number of parameters that should be optimized, i.e.,
         * the number of components of the function argument.
         *
         * @return number of components
         */
        virtual std::size_t dim() const = 0;

        /**
         * @brief General-purpose callback function
         *
         * This function is called after every iteration and can be used to
         * implement custom behavior, e.g., writing information to a log file,
         * producing visualizations, or collecting statistics. If extrapolation
         * is used, then it is also called after each successful extrapolation
         * step, so twice per iteration. A Boolean flag is used to distinguish
         * these two situations.
         *
         * @param iteration     current number of performed iterations
         * @param point         current vector of parameters
         * @param value         current function value
         * @param gradient      function gradient at current position
         * @param extrapolation true after extrapolation steps
         */
        virtual void hook(std::size_t iteration, const Point& point, Real value,
            const Point& gradient, bool extrapolation = false) const = 0;

        virtual ~ProblemBase(){};
      };

    /**
     * @brief Abstract base class for numerical differentiation
     *
     * Derive from this class to define an optimization problem that
     * uses numerical derivative evaluation based on finite differences.
     * Note that this doesn't scale very well to higher dimensions: in
     * \f$n\f$ dimensions, \f$n+1\f$ function evaluations are needed per
     * gradient evaluation. Try to supply the gradient as a closed
     * formulation, through automatic differentiation, or through an
     * adjoint problem instead, if possible.
     *
     * @tparam Real  data type of scalars and function values
     * @tparam Point data type of parameter vectors
     */
    template<typename Real, typename Point>
      class FiniteDifferenceProblemBase
      : public ProblemBase<Real,Point>
      {
        std::vector<Real> eps;
        mutable Point eval_point;
        mutable Real dir_deriv_val = 0.;

        public:

        /**
         * @brief Constructor for uniform step width
         *
         * @param eps_ step width for numerical gradient evaluation
         */
        FiniteDifferenceProblemBase(Real eps_ = 1e-10)
          : eps(1,eps_)
        {
          if (eps[0] <= 0.)
            throw std::logic_error("step width must be positive");
        }

        /**
         * @brief Constructor for different step widths per dimension
         *
         * @param eps_ step widths for numerical gradient evaluation
         */
        explicit FiniteDifferenceProblemBase(const Point& eps_)
        {
          for (std::size_t i = 0; i < eps_.size(); i++)
          {
            if (eps_[i] <= 0.)
              throw std::logic_error("step width must be positive");
            eps.push_back(eps_[i]);
          }
        }

        Real dir_deriv(const Point& point, const Point& direction,
            Point& gradient, bool subsequent = false) const override
        {
          const Real dir_eps = this->dir_eps(direction);

          eval_point = point;
          eval_point.axpy(direction,dir_eps);
          const Real plus_val = this->value(eval_point);

          dir_deriv_val = this->value(point);

          return (plus_val - dir_deriv_val) / dir_eps;
        }

        void gradient(const Point& point, Point& gradient,
            const Point* const direction = nullptr,
            Real derivative = 0.) const override
        {
          std::size_t skip = point.size();
          if (direction)
            skip = abs_max(*direction);

          eval_point = point;
          if (gradient.size() != this->dim())
            gradient = this->zero();

          for (std::size_t i = 0; i < point.size(); i++)
            if (i != skip)
            {
              const Real dir_eps = eps[i % eps.size()];
              eval_point[i] += dir_eps;
              const Real plus_val = this->value(eval_point);
              eval_point[i] -= dir_eps;

              gradient[i] = plus_val / dir_eps;
            }

          const Real val = direction ? dir_deriv_val : this->value(point);

          for (std::size_t i = 0; i < point.size(); i++)
            if (i != skip)
            {
              const Real dir_eps = eps[i % eps.size()];
              gradient[i] -= val / dir_eps;
            }

          if (skip < point.size())
          {
            gradient[skip] = derivative;

            for (std::size_t i = 0; i < point.size(); i++)
              if (i != skip)
                gradient[skip] -= gradient[i] * (*direction)[i];

            gradient[skip] /= (*direction)[skip];
          }
        }

        std::pair<Real,Real> value_and_deriv(const Point& point,
            const Point& direction, Point& gradient) const override
        {
          const Real dir_deriv = this->dir_deriv(point,direction,gradient);

          return {dir_deriv_val, dir_deriv};
        }

        private:

        /**
         * @brief Compute step width for directional derivatives
         *
         * If different step widths for numerical differentiation
         * per dimension are used, then a consistent choice has to
         * be made along search directions. This function works by
         * spanning a hyperellipsoid of step widths to obtain a
         * weighted average for the given search direction.
         *
         * @param direction direction for derivative computation
         *
         * @return step width in given direction
         */
        Real dir_eps(const Point& direction) const
        {
          if (eps.size() == 1)
            return eps[0] / std::sqrt(direction * direction);

          Real sum = 0.;
          for (std::size_t i = 0; i < direction.size(); i++)
            sum += eps[i]*eps[i] * direction[i]*direction[i];
          return std::sqrt(sum) / (direction * direction);
        }

        /**
         * @brief Determine component with maximum absolute value
         *
         * When computing the gradient using finite differences,
         * the previously computed directional derivative may be
         * reused, saving one function evaluation. This function
         * determines the dimension that has the largest contribution
         * to the search direction, which is the one that is skipped.
         *
         * @param direction direction used in derivative computation
         *
         * @return index of dimension that should be skipped
         */
        std::size_t abs_max(const Point& direction) const
        {
          std::size_t max_index = 0;
          Real max_val = 0.;
          for (unsigned int i = 0; i < direction.size(); i++)
          {
            const Real abs_val = std::abs(direction[i]);
            if (abs_val > max_val)
            {
              max_val = abs_val;
              max_index = i;
            }
          }

          return max_index;
        }
      };

#if HAVE_MUPARSER
    /**
     * @brief Defines an optimization problem by parsing an expression
     *
     * This class provides a simple way to define objective functions
     * that should be optimized. The function is passed as an
     * expression contained in a string, which is then parsed by
     * muParser. Up to 128 variables, and therefore dimensions, may be
     * introduced. The function gradient is evaluated through numerical
     * differences, with configurable step width.
     */
    class ParsedProblem
      : public FiniteDifferenceProblemBase<double,VectorClass<double>>
    {
      public:

        using Real = double;
        using Point = VectorClass<double>;

      private:

        mu::Parser parser;
        std::map<std::string,Real*> var;
        mutable Real val;
        mutable unsigned int eval_count = 0;

      public:

        /**
         * @brief Constructor for uniform step width
         *
         * @param expression string containing function definition
         * @param eps        step width for numerical gradient evaluation
         */
        ParsedProblem(const std::string& expression, Real eps = 1e-10)
          : FiniteDifferenceProblemBase(eps)
        {
          parser.SetVarFactory(add_variable,nullptr);
          parser.SetExpr(expression);
          std::cout << "expression: " << expression << std::endl;
          var = parser.GetUsedVar();
        }

        /**
         * @brief Constructor for different step widths per dimension
         *
         * @param expression string containing function definition
         * @param eps        step widths for numerical gradient evaluation
         */
        ParsedProblem(const std::string& expression, const Point& eps)
          : FiniteDifferenceProblemBase(eps)
        {
          parser.SetVarFactory(add_variable,nullptr);
          parser.SetExpr(expression);
          std::cout << "expression: " << expression << std::endl;
          var = parser.GetUsedVar();
        }

        Real value(const Point& point, bool subsequent = false) const override
        {
          if (subsequent)
            return val;

          if (point.size() != var.size())
            throw std::domain_error("x and var have different dimensions");

          std::size_t count = 0;
          for (auto& entry : var)
            *(entry.second) = point[count++];

          try
          {
            val = parser.Eval();
          }
          catch(mu::Parser::exception_type& e)
          {
            std::cout << e.GetMsg() << std::endl;
            val = std::numeric_limits<Real>::infinity();
          }

          eval_count++;
          return val;
        }

        Point zero() const override
        { 
          return Point(var.size(),0.);
        }

        std::size_t dim() const override
        {
          return var.size();
        }

        void hook(std::size_t iteration, const Point& point, Real value,
            const Point& gradient, bool extrapolation = false) const override
        {}

        unsigned int evals() const
        {
          return eval_count;
        }

      private:

        /**
         * @brief Function that provides storage space for muParser
         *
         * To pass function arguments to muParser, the storage location
         * has to be made known, which is done using this function.
         * Afterwards, coordinates may simply be written to these locations
         * calling the parser during function evaluations.
         *
         * @param str name of variable that requires storage
         * @param ptr custom user data (not used)
         *
         * @return address where the variable is stored
         */
        static double* add_variable(const char* str, void* ptr)
        {
          static Real var[128];
          static int index = -1;

          index++;

          std::cout << "Generating new variable \"" 
            << str << "\" (slots left: " 
            << 127-index << ")" << std::endl;

          var[index] = 0.;

          if (index > 127)
            throw mu::ParserError("Variable buffer overflow.");
          else
            return &var[index];
        }
    };
#endif // HAVE_MUPARSER

    /**
     * @brief Abstract base class for preconditioners
     *
     * Derive from this class to provide a preconditioner. The derived
     * class can be handed to the optimization method, which will then
     * apply the preconditioner in each iteration.
     */
    template<typename Real, typename Point>
      struct PreconditionerBase
      {
        /**
         * @brief Unique name to identify preconditioner
         *
         * @return string containing the name
         */
        virtual std::string name() const = 0;

        /**
         * @brief Apply preconditioner to gradient
         *
         * This method represents the application of the preconditioner
         * to a given vector. A copy of the current gradient vector is
         * handed to the preconditioner class, and it should modify its
         * function argument to return the result of applying the
         * preconditioner to the gradient.
         *
         * @param         point     current position
         * @param[in,out] prec_grad gradient, becomes preconditioned gradient
         */
        virtual void apply(const Point& point, Point& prec_grad) const = 0;

        virtual ~PreconditionerBase(){};
      };

  }
}

#endif // DUNE_NONLINOPT_USERINTERFACE_HH
