#pragma once

#include <optional>
#include <vector>

#include <dune/common/dynvector.hh>
#include <dune/common/version.hh>
#include <dune/localfunctions/lagrange.hh>
#include <dune/vtk/gridfunctions/common.hh>
#include <dune/vtk/gridfunctions/continuousgridfunction.hh>
#include <dune/vtk/utility/errors.hh>
#include <dune/vtk/utility/lagrangepoints.hh>

namespace Dune
{
  namespace Vtk
  {
    template <class GridType>
    class LagrangeGridCreator;

    /// \brief Grid-function representing values from a VTK file with local Lagrange
    /// interpolation of the values stored on the Lagrange nodes.
    template <class GridType, class FieldType, class Context>
    class LagrangeGridFunction
    {
      using Grid = GridType;
      using Field = FieldType;

      using Factory = GridFactory<Grid>;
      using GridCreator = LagrangeGridCreator<Grid>;

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
        using LFE = LagrangeLocalFiniteElement<LagrangePointSet, LC::dimension, FieldType, FieldType>;
        using LB = typename LFE::Traits::LocalBasisType;

      public:
        using LocalContext = LC;
        using Domain = typename LC::Geometry::LocalCoordinate;
        using Range = DynamicVector<Field>;
        using Signature = Range(Domain);

      public:
        /// Constructor. Stores references to the passed data.
        PointDataLocalFunction (GridCreator const* creator, std::vector<Field> const* values, unsigned int comp,
                                std::vector<std::uint8_t> const* types,
                                std::vector<std::int64_t> const* offsets,
                                std::vector<std::int64_t> const* connectivity)
          : creator_(creator)
          , values_(values)
          , comp_(comp)
          , types_(types)
          , offsets_(offsets)
          , connectivity_(connectivity)
        {}

        PointDataLocalFunction () = default;

        /// Binding the local-function to an element.
        /**
         * Constructs a new local finite-element with a polynomial order given
         * by the number of Lagrange nodes on the element. Extracts values on all
         * Lagrange nodes stored in file and stores these local coefficients in
         * a class member variable.
         **/
        void bind (LocalContext const& element)
        {
          unsigned int insertionIndex = creator_->factory().insertionIndex(element);

          std::int64_t shift = (insertionIndex == 0 ? 0 : (*offsets_)[insertionIndex-1]);
          std::int64_t numNodes = (*offsets_)[insertionIndex] - shift;
          [[maybe_unused]] std::int64_t maxNumNodes = numLagrangePoints(element.type(), 20);
          VTK_ASSERT(numNodes > 0 && numNodes < maxNumNodes);

          int order = creator_->order(element.type(), numNodes);
          VTK_ASSERT(order > 0 && order < 20);

          // construct a local finite-element with the corresponding order and Lagrange points
          // as stored in the file
          lfe_.emplace(LFE{element.type(), (unsigned int)(order)});
          lgeo_.emplace(creator_->localGeometry(element));

          // collect values on lagrange nodes
          localValues_.resize(numNodes);
          for (std::int64_t i = shift, i0 = 0; i < (*offsets_)[insertionIndex]; ++i, ++i0) {
            std::int64_t idx = (*connectivity_)[i];
            DynamicVector<Field>& v = localValues_[i0];
            v.resize(comp_);
            for (unsigned int j = 0; j < comp_; ++j)
              v[j] = (*values_)[comp_*idx + j];
          }
        }

        /// Unbind the local-function and the local finite-element from the element.
        void unbind ()
        {
          lfe_.reset();
          lgeo_.reset();
        }

        /// Evaluation in element local coordinates. Essentially, a local Lagrange
        /// interpolation with coefficients extracted in \ref bind().
        // NOTE: do we need to transform the local coordinates?
        Range operator() (Domain const& local) const
        {
          assert(!!lfe_);
          auto const& lb = lfe_->localBasis();
          lb.evaluateFunction(lgeo_->global(local), shapeValues_);
          assert(shapeValues_.size() == localValues_.size());

          Range y(comp_, Field(0));
          for (std::size_t i = 0; i < shapeValues_.size(); ++i)
            y.axpy(shapeValues_[i], localValues_[i]);

          return y;
        }

      private:
        GridCreator const* creator_ = nullptr;
        std::vector<Field> const* values_ = nullptr;
        unsigned int comp_;
        std::vector<std::uint8_t> const* types_ = nullptr;
        std::vector<std::int64_t> const* offsets_ = nullptr;
        std::vector<std::int64_t> const* connectivity_ = nullptr;

        // Local Finite-Element
        std::optional<LFE> lfe_ = std::nullopt;
        std::optional<typename GridCreator::LocalGeometry> lgeo_ = std::nullopt;

        // cache of local values
        std::vector<DynamicVector<Field>> localValues_;
        mutable std::vector<typename LB::Traits::RangeType> shapeValues_;
      };

      /// Evaluation of data on the cells of the grid
      template <class LC>
      class CellDataLocalFunction
      {
      public:
        using LocalContext = LC;
        using Domain = typename LC::Geometry::LocalCoordinate;
        using Range = DynamicVector<Field>;
        using Signature = Range(Domain);

      public:
        /// Constructor. Stores references to the passed data.
        CellDataLocalFunction (GridCreator const* creator, std::vector<Field> const* values, unsigned int comp,
                                std::vector<std::uint8_t> const* /*types*/,
                                std::vector<std::int64_t> const* /*offsets*/,
                                std::vector<std::int64_t> const* /*connectivity*/)
          : creator_(creator)
          , values_(values)
          , comp_(comp)
        {}

        CellDataLocalFunction () = default;

        /// Binding the local-function to an element extract the cell-value from the vector
        /// of data.
        void bind (LocalContext const& element)
        {
          unsigned int idx = creator_->factory().insertionIndex(element);

          // collect values on cells
          DynamicVector<Field>& v = localValue_;
          v.resize(comp_);

          for (unsigned int j = 0; j < comp_; ++j)
            v[j] = (*values_)[comp_*idx + j];
        }

        /// Unbinds from the bound element. Does nothing
        void unbind ()
        {}

        /// Evaluation in local element coordinates. Returns the constant value
        /// extracted in \ref bind().
        Range operator() (Domain const& local) const
        {
          return localValue_;
        }

      private:
        GridCreator const* creator_ = nullptr;
        std::vector<Field> const* values_ = nullptr;
        unsigned int comp_;

        // cache of local values
        DynamicVector<Field> localValue_;
      };

      /// Type switch for local-functions depending on the Context
      template <class LC>
      using LocalFunction = std::conditional_t< std::is_same_v<Context,PointContext>,
        PointDataLocalFunction<LC>,
        CellDataLocalFunction<LC>>;

    public:
      /// Construct a grid-function. Passed in data is stroed by reference, thus must have
      /// a life-time greater than that of the grid-function and corresponding local-function.
      LagrangeGridFunction (GridCreator const& creator, std::vector<Field> const& values,
                            std::string name ,unsigned int ncomps, Vtk::DataTypes dataType,
                            std::vector<std::uint8_t> const& types,
                            std::vector<std::int64_t> const& offsets,
                            std::vector<std::int64_t> const& connectivity)
        : creator_(&creator)
        , values_(&values)
        , name_(std::move(name))
        , ncomps_(ncomps)
        , dataType_(dataType)
        , types_(&types)
        , offsets_(&offsets)
        , connectivity_(&connectivity)
      {}

      LagrangeGridFunction () = default;

      /// Global evaluation. Not supported!
      Range operator() (Domain const& global) const
      {
        DUNE_THROW(Dune::NotImplemented, "Evaluation in global coordinates not implemented.");
        return Range(ncomps_, 0);
      }

      /// Return a type that defines the element that can be iterated.
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

      /// Construct a local-function depending on the Context type either PointDataLocalFunction
      /// or CellDataLocalFunction
      friend LocalFunction<typename EntitySet::Element> localFunction (LagrangeGridFunction const& gf)
      {
        return {gf.creator_, gf.values_, gf.ncomps_, gf.types_, gf.offsets_, gf.connectivity_};
      }

    private:
      GridCreator const* creator_ = nullptr;
      std::vector<Field> const* values_ = nullptr;
      std::string name_ = "GridFunction";
      unsigned int ncomps_ = 0;
      Vtk::DataTypes dataType_ = Vtk::DataTypes::UNKNOWN;
      std::vector<std::uint8_t> const* types_ = nullptr;
      std::vector<std::int64_t> const* offsets_ = nullptr;
      std::vector<std::int64_t> const* connectivity_ = nullptr;

      EntitySet entitySet_;
    };

  } // end namespace Vtk
} // end namespace Dune
