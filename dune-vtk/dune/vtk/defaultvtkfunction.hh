#pragma once

#include <type_traits>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

#include "localfunctioninterface.hh"

namespace Dune
{
  namespace Vtk
  {
    /// Type erasure for dune-functions LocalFunction interface
    template <class GridView, class LocalFunction>
    class LocalFunctionWrapper final
        : public LocalFunctionInterface<GridView>
    {
      using Self = LocalFunctionWrapper;
      using Interface = LocalFunctionInterface<GridView>;
      using Entity = typename Interface::Entity;
      using LocalCoordinate = typename Interface::LocalCoordinate;

    public:
      /// Constructor. Stores a copy of the passed `localFct` in a local variable.
      template <class LocalFct,
        disableCopyMove<Self, LocalFct> = 0>
      explicit LocalFunctionWrapper (LocalFct&& localFct)
        : localFct_(std::forward<LocalFct>(localFct))
      {}

      /// Bind the LocalFunction to the Entity
      virtual void bind (Entity const& entity) override
      {
        localFct_.bind(entity);
      }

      /// Unbind the LocalFunction from the Entity
      virtual void unbind () override
      {
        localFct_.unbind();
      }

      /// Evaluate the LocalFunction in LocalCoordinates
      virtual double evaluate (int comp, LocalCoordinate const& xi) const override
      {
        return evaluateImpl(comp, localFct_(xi));
      }

    private:
      // Evaluate a component of a vector valued data
      template <class T, int N, int M>
      double evaluateImpl (int comp, FieldMatrix<T,N,M> const& mat) const
      {
        int r = comp / 3;
        int c = comp % 3;
        return r < N && c < M ? mat[r][c] : 0.0;
      }

      // Evaluate a component of a vector valued data
      template <class T, int N>
      double evaluateImpl (int comp, FieldVector<T,N> const& vec) const
      {
        return comp < N ? vec[comp] : 0.0;
      }

      // Evaluate a component of a vector valued data
      template <class T,
        std::enable_if_t<IsIndexable<T,int>::value, int> = 0>
      double evaluateImpl (int comp, T const& value) const
      {
        return value[comp];
      }

      // Return the scalar values
      template <class T,
        std::enable_if_t<not IsIndexable<T,int>::value, int> = 0>
      double evaluateImpl (int comp, T const& value) const
      {
        assert(comp == 0);
        return value;
      }

    private:
      LocalFunction localFct_;
    };

  } // end namespace Vtk
} // end namespace Dune
