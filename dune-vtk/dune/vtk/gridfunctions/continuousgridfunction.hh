#pragma once

#include <vector>

#include <dune/common/dynvector.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>
#include <dune/vtk/gridfunctions/common.hh>

namespace Dune
{
  namespace Vtk
  {
    /// \brief A GridFunction representing data stored on the grid vertices in a continuous manner.
    template <class GridType, class FieldType, class Context>
    class ContinuousGridFunction
    {
      using Grid = GridType;
      using Field = FieldType;

      using Factory = GridFactory<Grid>;

    public:
      struct EntitySet
      {
        using Grid = GridType;
        using Element = typename GridType::template Codim<0>::Entity;
        using LocalCoordinate = typename Element::Geometry::LocalCoordinate;
        using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;
      };

      using Domain = typename EntitySet::GlobalCoordinate;
      using Range = DynamicVector<Field>;
      using Signature = Range(Domain);

    private:
      template <class LC>
      class PointDataLocalFunction
      {
        // NOTE: local finite-element fixed to Lagrange P1/Q1
        using LFECache = LagrangeLocalFiniteElementCache<Field,Field,LC::mydimension,1>;
        using LFE = typename LFECache::FiniteElementType;
        using LB = typename LFE::Traits::LocalBasisType;

      public:
        using LocalContext = LC;
        using Domain = typename LC::Geometry::LocalCoordinate;
        using Range = DynamicVector<Field>;
        using Signature = Range(Domain);

      public:
        PointDataLocalFunction (Factory const* factory, std::vector<Field> const* values, unsigned int comp)
          : factory_(factory)
          , values_(values)
          , comp_(comp)
        {}

        PointDataLocalFunction () = default;

        void bind (LocalContext const& element)
        {
          lfe_ = &cache_.get(element.type());

          // collect values on vertices
          // NOTE: assumes, that Lagrange nodes are ordered like element vertices
          localValues_.resize(element.subEntities(Grid::dimension));
          for (unsigned int i = 0; i < element.subEntities(Grid::dimension); ++i) {
            unsigned int idx = factory_->insertionIndex(element.template subEntity<Grid::dimension>(i));
            DynamicVector<Field>& v = localValues_[i];
            v.resize(comp_);
            for (unsigned int j = 0; j < comp_; ++j)
              v[j] = (*values_)[comp_*idx + j];
          }
        }

        void unbind ()
        {
          lfe_ = nullptr;
        }

        Range operator() (Domain const& local) const
        {
          assert(!!lfe_);
          auto const& lb = lfe_->localBasis();
          lb.evaluateFunction(local, shapeValues_);
          assert(shapeValues_.size() == localValues_.size());

          Range y(comp_, Field(0));
          for (std::size_t i = 0; i < shapeValues_.size(); ++i)
            y.axpy(shapeValues_[i], localValues_[i]);

          return y;
        }

      private:
        Factory const* factory_ = nullptr;
        std::vector<Field> const* values_ = nullptr;
        unsigned int comp_;

        // Local Finite-Element
        LFECache cache_;
        LFE const* lfe_ = nullptr;

        // cache of local values
        std::vector<DynamicVector<Field>> localValues_;
        mutable std::vector<typename LB::Traits::RangeType> shapeValues_;
      };

      template <class LC>
      class CellDataLocalFunction
      {
      public:
        using LocalContext = LC;
        using Domain = typename LC::Geometry::LocalCoordinate;
        using Range = DynamicVector<Field>;
        using Signature = Range(Domain);

      public:
        CellDataLocalFunction (Factory const* factory, std::vector<Field> const* values, unsigned int comp)
          : factory_(factory)
          , values_(values)
          , comp_(comp)
        {}

        CellDataLocalFunction () = default;

        void bind (LocalContext const& element)
        {
          unsigned int idx = factory_->insertionIndex(element);

          // collect values on cells
          DynamicVector<Field>& v = localValue_;
          v.resize(comp_);

          for (unsigned int j = 0; j < comp_; ++j)
            v[j] = (*values_)[comp_*idx + j];
        }

        void unbind ()
        {}

        Range operator() (Domain const& local) const
        {
          return localValue_;
        }

      private:
        Factory const* factory_ = nullptr;
        std::vector<Field> const* values_ = nullptr;
        unsigned int comp_;

        // cache of local values
        DynamicVector<Field> localValue_;
      };

      template <class LC>
      using LocalFunction = std::conditional_t< std::is_same_v<Context,PointContext>,
        PointDataLocalFunction<LC>,
        CellDataLocalFunction<LC>>;

    public:
      template <class GridCreator>
      ContinuousGridFunction (GridCreator const& creator, std::vector<Field> const& values,
                              std::string name, unsigned int ncomps, Vtk::DataTypes dataType,
                              std::vector<std::uint8_t> const& /*types*/,
                              std::vector<std::int64_t> const& /*offsets*/,
                              std::vector<std::int64_t> const& /*connectivity*/)
        : factory_(&creator.factory())
        , values_(&values)
        , name_(std::move(name))
        , ncomps_(ncomps)
        , dataType_(dataType)
      {}

      ContinuousGridFunction () = default;

      Range operator() (Domain const& global) const
      {
        DUNE_THROW(Dune::NotImplemented, "Evaluation in global coordinates not implemented.");
        return Range(ncomps_, 0);
      }

      EntitySet const& entitySet () const
      {
        return entitySet_;
      }

      std::string const& name () const
      {
        return name_;
      }

      int numComponents () const
      {
        return ncomps_;
      }

      Vtk::DataTypes dataType () const
      {
        return dataType_;
      }

      friend LocalFunction<typename EntitySet::Element> localFunction (ContinuousGridFunction const& gf)
      {
        return {gf.factory_, gf.values_, gf.ncomps_};
      }

    private:
      Factory const* factory_;
      std::vector<Field> const* values_ = nullptr;
      std::string name_ = "GridFunction";
      unsigned int ncomps_ = 0;
      Vtk::DataTypes dataType_ = Vtk::DataTypes::UNKNOWN;

      EntitySet entitySet_;
    };

  } // end namespace Vtk
} // end namespace Dune
