#pragma once

#include <memory>

#include <dune/grid/io/file/vtk/function.hh>

#include "localfunctioninterface.hh"

namespace Dune
{
  namespace Vtk
  {
    /// Type erasure for Legacy VTKFunction
    template <class GridView>
    class VTKLocalFunctionWrapper final
        : public LocalFunctionInterface<GridView>
    {
      using Interface = LocalFunctionInterface<GridView>;
      using Entity = typename Interface::Entity;
      using LocalCoordinate = typename Interface::LocalCoordinate;

    public:
      /// Constructor. Stores a shared pointer to the passed Dune::VTKFunction
      explicit VTKLocalFunctionWrapper (std::shared_ptr<VTKFunction<GridView> const> const& fct)
        : fct_(fct)
      {}

      /// Stores a pointer to the passed entity
      virtual void bind (Entity const& entity) override
      {
        entity_ = &entity;
      }

      /// Unsets the stored entity pointer
      virtual void unbind () override
      {
        entity_ = nullptr;
      }

      /// Evaluate the Dune::VTKFunction in LocalCoordinates on the stored Entity
      virtual double evaluate (int comp, LocalCoordinate const& xi) const override
      {
        return fct_->evaluate(comp, *entity_, xi);
      }

    private:
      std::shared_ptr<VTKFunction<GridView> const> fct_;
      Entity const* entity_;
    };

  } //end namespace Vtk
} // end namespace Dune
