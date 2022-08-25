#ifndef DUNE_NONLINOPT_VECTORCLASS_HH
#define DUNE_NONLINOPT_VECTORCLASS_HH

#include<array>
#include<vector>

#if HAVE_DUNE_COMMON
#include<dune/common/exceptions.hh>
#endif // HAVE_DUNE_COMMON

namespace Dune {
  namespace NonlinOpt {

    /**
     * @brief A simple numerical vector class
     *
     * This class can be used as the basis of optimization
     * problem definitions if there are no special requirements,
     * i.e., if the program is sequential and a flat array of
     * real numbers is an appropriate representation of the
     * parameters. Parallel computations are straightforward and
     * only require replacing this class with an alternative that
     * provides parallelized numerical linear algebra and scalar
     * products.
     */
    template<typename Real>
      class VectorClass
      {
        std::vector<Real> values;

        public:

        /**
         * @brief Default constructor
         *
         * Produces an empty vector that will later be resized
         * to fit the dimension of the optimization problem.
         */
        VectorClass()
        {}

        /**
         * @brief Fixed-size constructor
         *
         * Produces a vector of given size, optionally initialized
         * to a certain value in each component.
         *
         * @param n     number of components
         * @param value initial value of components
         */
        VectorClass(std::size_t n, Real value = 0.)
          : values(n,value)
        {}

        /**
         * @brief Constructor based on list of components
         *
         * Produces a vector from a given brace-enclosed list of
         * components.
         *
         * @param l list of vector components
         */
        VectorClass(const std::initializer_list<Real>& l)
          : values(l)
        {}

        /**
         * @brief Scaled addition method
         *
         * Adds a multiple of another vector to the current one,
         * i.e., \f$y = a * x + y\f$.
         *
         * @param other vector that will be scaled and added
         * @param alpha scaling factor
         */
        void axpy(const VectorClass& other, Real alpha)
        {
          if (values.size() != other.values.size())
          {
#if HAVE_DUNE_COMMON
            DUNE_THROW(Dune::Exception,"size mismatch");
#else
            throw std::logic_error("size mismatch");
#endif // HAVE_DUNE_COMMON
          }

          const std::size_t n = values.size();
          for (std::size_t i = 0; i < n; i++)
            values[i] += alpha * other.values[i];
        }

        /**
         * @brief Vector scaling
         *
         * Multiplies each component with a given scalar.
         *
         * @param alpha scaling factor
         */
        void operator*=(Real alpha)
        {
          for (std::size_t i = 0; i < values.size(); i++)
            values[i] *= alpha;
        }

        /**
         * @brief Vector addition
         *
         * Adds another vector to the current one.
         *
         * @param other vector to add
         */
        void operator+=(const VectorClass& other)
        {
          if (values.size() != other.values.size())
          {
#if HAVE_DUNE_COMMON
            DUNE_THROW(Dune::Exception,"size mismatch");
#else
            throw std::logic_error("size mismatch");
#endif // HAVE_DUNE_COMMON
          }

          for (std::size_t i = 0; i < values.size(); i++)
            values[i] += other.values[i];
        }

        /**
         * @brief Vector subtraction
         *
         * Subtracts another vector from the current one.
         *
         * @param other vector to subtract
         */
        void operator-=(const VectorClass& other)
        {
          if (values.size() != other.values.size())
          {
#if HAVE_DUNE_COMMON
            DUNE_THROW(Dune::Exception,"size mismatch");
#else
            throw std::logic_error("size_mismatch");
#endif // HAVE_DUNE_COMMON
          }

          for (std::size_t i = 0; i < values.size(); i++)
            values[i] -= other.values[i];
        }

        /**
         * @brief Dimension of the vector
         *
         * @return number of components
         */
        std::size_t size() const
        {
          return values.size();
        }

        /**
         * @brief Read-only element access
         *
         * @param i index of desired component
         *
         * @return value of i-th component
         */
        Real operator[](std::size_t i) const
        {
          return values[i];
        }

        /**
         * @brief Read-write element access
         *
         * @param i index of desired component
         *
         * @return reference to i-th component
         */
        Real& operator[](std::size_t i)
        {
          return values[i];
        }

        /**
         * @brief Scalar product
         *
         * Produces scalar product with other given vector.
         *
         * @param other vector to multiply the current one with
         *
         * @return sum of component-wise products
         */
        Real operator*(const VectorClass& other) const
        {
          std::array<Real,4> out = {0.,0.,0.,0.};

          const std::size_t n = values.size();
          const std::size_t n4 = n % 4;
          for (std::size_t i = 0; i < n4; i++)
            out[i] += values[i] * other.values[i];
          for (std::size_t i = n4; i < n; i += 4)
          {
            out[0] += values[i]   * other.values[i];
            out[1] += values[i+1] * other.values[i+1];
            out[2] += values[i+2] * other.values[i+2];
            out[3] += values[i+3] * other.values[i+3];
          }

          return out[0]+out[1]+out[2]+out[3];
        }

        /**
         * @brief Maximum norm
         *
         * @return maximum absolute value of components
         */
        Real inf_norm() const
        {
          Real max_value = 0.;
          for (std::size_t i = 0; i < values.size(); i++)
          {
            const Real abs_value = std::abs(values[i]);
            if (abs_value > max_value)
              max_value = abs_value;
          }

          return max_value;
        }
      };

  }
}

#endif // DUNE_NONLINOPT_VECTORCLASS_HH
