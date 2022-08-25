#pragma once

#include <memory>
#include <type_traits>

#include <dune/common/typetraits.hh>

#include "localfunctioninterface.hh"
#include "legacyvtkfunction.hh"
#include "defaultvtkfunction.hh"

namespace Dune
{
  namespace Vtk
  {
    /// \brief A Vtk::LocalFunction is a function-like object that can be bound to a grid element
    /// an that provides an evaluate method with a component argument.
    /**
    * Stores internally a Vtk::LocalFunctionInterface object for the concrete evaluation.
    **/
    template <class GridView>
    class LocalFunction
    {
      using Self = LocalFunction;
      using Entity = typename GridView::template Codim<0>::Entity;
      using LocalCoordinate = typename Entity::Geometry::LocalCoordinate;

      template <class LF, class E>
      using HasBind = decltype((std::declval<LF>().bind(std::declval<E>()), true));

    private:
      struct RangeProxy
      {
        using value_type = double;
        using field_type = double;

        RangeProxy (LocalFunctionInterface<GridView> const& localFct,
                    std::vector<int> const& components,
                    LocalCoordinate const& local)
          : localFct_(localFct)
          , components_(components)
          , local_(local)
        {}

        std::size_t size () const
        {
          return components_.size();
        }

        double operator[] (std::size_t i) const
        {
          return i < size() ? localFct_.evaluate(components_[i], local_) : 0.0;
        }

      private:
        LocalFunctionInterface<GridView> const& localFct_;
        std::vector<int> const& components_;
        LocalCoordinate local_;
      };

    public:
      /// Construct the Vtk::LocalFunction from any function object that has a bind(element) method.
      template <class LF,
        disableCopyMove<Self, LF> = 0,
        HasBind<LF,Entity> = true>
      explicit LocalFunction (LF&& lf)
        : localFct_(std::make_shared<LocalFunctionWrapper<GridView,LF>>(std::forward<LF>(lf)))
      {}

      /// Construct a Vtk::LocalFunction from a legacy VTKFunction
      explicit LocalFunction (std::shared_ptr<VTKFunction<GridView> const> const& lf)
        : localFct_(std::make_shared<VTKLocalFunctionWrapper<GridView>>(lf))
      {}

      /// Allow the default construction of a Vtk::LocalFunction. After construction, the
      /// LocalFunction is in an invalid state.
      LocalFunction () = default;

      /// Bind the function to the grid entity
      void bind (Entity const& entity)
      {
        assert(bool(localFct_));
        localFct_->bind(entity);
      }

      /// Unbind from the currently bound entity
      void unbind ()
      {
        assert(bool(localFct_));
        localFct_->unbind();
      }

      /// Return a proxy object to access the components of the range vector
      RangeProxy operator() (LocalCoordinate const& xi) const
      {
        assert(bool(localFct_));
        return {*localFct_, components_, xi};
      }

      /// Evaluate the `c`th component of the Range value at local coordinate `xi`
      double evaluate (int c, LocalCoordinate const& xi) const
      {
        assert(bool(localFct_));
        return c < int(components_.size()) ? localFct_->evaluate(components_[c], xi) : 0.0;
      }

      void setComponents (std::vector<int> components)
      {
        components_ = std::move(components);
      }

    private:
      std::shared_ptr<LocalFunctionInterface<GridView>> localFct_ = nullptr;
      std::vector<int> components_;
    };

  } // end namespace Vtk
} // end namespace Dune
