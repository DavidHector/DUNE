#pragma once

namespace Dune
{
  namespace Vtk
  {
    /// \brief An abstract base class for LocalFunctions that can be bound to an element and
    /// evaluated in local coordinates w.r.t. to a component of its value.
    template <class GridView>
    class LocalFunctionInterface
    {
    public:
      using Entity = typename GridView::template Codim<0>::Entity;
      using LocalCoordinate = typename Entity::Geometry::LocalCoordinate;

      /// Bind the function to the grid entity
      virtual void bind (Entity const& entity) = 0;

      /// Unbind from the currently bound entity
      virtual void unbind () = 0;

      /// Evaluate single component comp in the entity at local coordinates xi
      virtual double evaluate (int comp, LocalCoordinate const& xi) const = 0;

      /// Virtual destructor
      virtual ~LocalFunctionInterface () = default;
    };

  } // end namespace Vtk
} // end namespace Dune
