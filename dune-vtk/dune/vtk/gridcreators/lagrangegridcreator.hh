#pragma once

#include <cassert>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/geometry/utility/typefromvertexcount.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/localfunctions/lagrange.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/vtk/types.hh>
#include <dune/vtk/gridcreatorinterface.hh>
#include <dune/vtk/gridfunctions/lagrangegridfunction.hh>
#include <dune/vtk/utility/lagrangepoints.hh>

namespace Dune
{
  namespace Vtk
  {
    // \brief Create a grid from data that represents higher (lagrange) cells.
    /**
    * The grid is created from the first nodes of a cell parametrization, representing
    * the  corner vertices. Thus a piecewise "flat" grid is constructed. The
    * parametrization is 1. passed as a local element parametrization to the
    * `insertElement()` function of a gridFactory to allow the grid itself to handle the
    * parametrization and 2. is stored internally that can be accessed by using this
    * GridCreator object as a grid function, or by extracting locally the parametrization
    * on each existing grid element after creation of the grid.
    *
    * So, the LagrangeGridCreator models both, a `GridCreator` and a `GridFunction`.
    **/
    template <class GridType>
    struct LagrangeGridCreator
        : public GridCreatorInterface<GridType, LagrangeGridCreator<GridType>>
    {
      using Self = LagrangeGridCreator;
      using Super = GridCreatorInterface<GridType, Self>;
      using GlobalCoordinate = typename Super::GlobalCoordinate;

      using Nodes = std::vector<GlobalCoordinate>;

      struct ElementParametrization
      {
        GeometryType type;                  //< Geometry type of the element
        std::vector<std::int64_t> nodes;    //< Indices of the w.r.t. `nodes_` vector
        std::vector<unsigned int> corners;  //< Insertion-indices of the element corner nodes
      };

      using Parametrization = std::vector<ElementParametrization>;
      using Element = typename GridType::template Codim<0>::Entity;
      using LocalCoordinate = typename Element::Geometry::LocalCoordinate;

      class LocalParametrization;
      class LocalFunction;

    public:
      using LocalGeometry = MultiLinearGeometry<typename Element::Geometry::ctype,Element::dimension,Element::dimension>;

    public:
      using Super::Super;
      using Super::factory;

      /// Implementation of the interface function `insertVertices()`
      void insertVerticesImpl (std::vector<GlobalCoordinate> const& points,
                              std::vector<std::uint64_t> const& /*point_ids*/)
      {
        // store point coordinates in member variable
        nodes_ = points;
      }

      template <class F>
      using HasParametrizedElements = decltype(std::declval<F>().insertElement(std::declval<GeometryType>(),
        std::declval<std::vector<unsigned int> const&>(), std::declval<std::function<GlobalCoordinate(LocalCoordinate)>>()));

      /// Implementation of the interface function `insertElements()`
      void insertElementsImpl (std::vector<std::uint8_t> const& types,
                              std::vector<std::int64_t> const& offsets,
                              std::vector<std::int64_t> const& connectivity)
      {
        assert(nodes_.size() > 0);

        // mapping of node index to element-vertex index
        std::vector<std::int64_t> elementVertices(nodes_.size(), -1);
        parametrization_.reserve(types.size());

        std::int64_t vertexIndex = 0;
        for (std::size_t i = 0; i < types.size(); ++i) {
          auto type = Vtk::to_geometry(types[i]);
          if (type.dim() != GridType::dimension)
            continue;

          Vtk::CellType cellType{type};
          auto refElem = referenceElement<double,GridType::dimension>(type);

          std::int64_t shift = (i == 0 ? 0 : offsets[i-1]);
          int nNodes = offsets[i] - shift;
          int nVertices = refElem.size(GridType::dimension);

          // insert vertices into grid and construct element vertices
          std::vector<unsigned int> element(nVertices);
          for (int j = 0; j < nVertices; ++j) {
            auto index = connectivity.at(shift + j);
            auto& vertex = elementVertices.at(index);
            if (vertex < 0) {
              factory().insertVertex(nodes_.at(index));
              vertex = vertexIndex++;
            }
            element[j] = vertex;
          }

          // permute element indices
          if (!cellType.noPermutation()) {
            // apply index permutation
            std::vector<unsigned int> cell(element.size());
            for (std::size_t j = 0; j < element.size(); ++j)
              cell[j] = element[cellType.permutation(j)];
            std::swap(element, cell);
          }

          // fill vector of element parametrizations
          parametrization_.push_back(ElementParametrization{type});
          auto& param = parametrization_.back();

          param.nodes.resize(nNodes);
          for (int j = 0; j < nNodes; ++j)
            param.nodes[j] = connectivity.at(shift + j);
          param.corners = element;

          // try to create element with parametrization
          if constexpr (Std::is_detected_v<HasParametrizedElements, GridFactory<GridType>>) {
            try {
              factory().insertElement(type, element,
                localParametrization(parametrization_.size()-1));
            } catch (Dune::GridError const& /* notImplemented */) {
              factory().insertElement(type, element);
            }
          } else {
            factory().insertElement(type, element);
          }
        }
      }

      /// \brief Construct an element parametrization
      /**
      * The returned LocalParametrization is a mapping `GlobalCoordinate(LocalCoordinate)`
      * where `LocalCoordinate is w.r.t. the local coordinate system in an element with
      * given `insertionIndex` (defined by the inserted corner vertices) and
      * `GlobalCoordinate` a world coordinate in the parametrized grid.
      **/
      LocalParametrization localParametrization (unsigned int insertionIndex) const
      {
        assert(!nodes_.empty() && !parametrization_.empty());
        auto const& localParam = parametrization_.at(insertionIndex);
        return LocalParametrization{nodes_, localParam, order(localParam)};
      }

      /// \brief Construct an element parametrization
      /**
      * The returned LocalParametrization is a mapping `GlobalCoordinate(LocalCoordinate)`
      * where `LocalCoordinate is w.r.t. the local coordinate system in the passed element
      * and `GlobalCoordinate` a world coordinate in the parametrized grid.
      **/
      LocalParametrization localParametrization (Element const& element) const
      {
        VTK_ASSERT(!nodes_.empty() && !parametrization_.empty());

        unsigned int insertionIndex = factory().insertionIndex(element);
        auto const& localParam = parametrization_.at(insertionIndex);
        VTK_ASSERT(element.type() == localParam.type);

        return {nodes_, localParam, order(localParam), localGeometry(element, localParam)};
      }

      /// \brief Construct a transformation of local element coordinates
      /**
      * An element might have a different local coordinate system than the coordinate
      * system used to defined the element parametrization. Thus coordinate transform
      * of the local parametrization is needed for element-local evaluations. This
      * local geometry transform is obtained by figuring out the permutation of corners
      * in the element corresponding to the inserted corner vertices.
      **/
      LocalGeometry localGeometry (Element const& element) const
      {
        VTK_ASSERT(!nodes_.empty() && !parametrization_.empty());

        unsigned int insertionIndex = factory().insertionIndex(element);
        auto const& localParam = parametrization_.at(insertionIndex);
        VTK_ASSERT(element.type() == localParam.type);

        return localGeometry(element, localParam);
      }

    private:
      // implementation details of localGeometry()
      LocalGeometry localGeometry (Element const& element, ElementParametrization const& localParam) const
      {
        // collect indices of vertices
        std::vector<unsigned int> indices(element.subEntities(GridType::dimension));
        for (unsigned int i = 0; i < element.subEntities(GridType::dimension); ++i)
          indices[i] = factory().insertionIndex(element.template subEntity<GridType::dimension>(i));

        // calculate permutation vector
        std::vector<unsigned int> permutation(indices.size());
        for (std::size_t i = 0; i < indices.size(); ++i) {
          auto it = std::find(localParam.corners.begin(), localParam.corners.end(), indices[i]);
          VTK_ASSERT(it != localParam.corners.end());
          permutation[i] = std::distance(localParam.corners.begin(), it);
        }

        auto refElem = referenceElement<typename Element::Geometry::ctype,Element::dimension>(localParam.type);
        std::vector<LocalCoordinate> corners(permutation.size());
        for (std::size_t i = 0; i < permutation.size(); ++i)
          corners[i] = refElem.position(permutation[i], Element::dimension);

        return {localParam.type, corners};
      }

    public:

      /// Determine lagrange order from number of points
      int order (GeometryType type, std::size_t nNodes) const
      {
        for (int o = 1; o <= int(nNodes); ++o)
          if (numLagrangePoints(type, o) == std::size_t(nNodes))
            return o;

        return 1;
      }

      int order (ElementParametrization const& localParam) const
      {
        return order(localParam.type, localParam.nodes.size());
      }

      /// Determine lagrange order from number of points from the first element parametrization
      int order () const
      {
        assert(!parametrization_.empty());
        auto const& localParam = parametrization_.front();
        return order(localParam);
      }

    public:
      /// \brief Local function representing the parametrization of the grid.
      /**
      * The returned object models Functions::Concept::LocalFunction
      * and can thus be bound to an element of the created grid and evaluated in
      * the local coordinates of the bound element.
      *
      * It is implemented in terms of the \ref LocalParametrization function
      * returned by the method \ref localParametrization(element). See comments
      * there for further details.
      *
      * Note, this methods requires the GridCreator to be based by
      * lvalue-reference. This is necessary, since we want to guarantee that all
      * internal storage is preserved while evaluating the local function.
      **/
      friend LocalFunction localFunction (LagrangeGridCreator& gridCreator)
      {
        return LocalFunction{gridCreator};
      }

      friend LocalFunction localFunction (LagrangeGridCreator const& gridCreator)
      {
        return LocalFunction{gridCreator};
      }

      friend LocalFunction localFunction (LagrangeGridCreator&& gridCreator)
      {
        DUNE_THROW(Dune::Exception, "Cannot pass temporary LagrangeGridCreator to localFunction(). Pass an lvalue-reference instead.");
        return LocalFunction{gridCreator};
      }

      std::string name () const
      {
        return "LagrangeParametrization";
      }

      int numComponents () const
      {
        return GlobalCoordinate::size();
      }

      Vtk::DataTypes dataType () const
      {
        return dataTypeOf<typename Element::Geometry::ctype>();
      }

      struct EntitySet
      {
        using Grid = GridType;
        using GlobalCoordinate = typename Self::GlobalCoordinate;
      };

      /// Dummy function returning a placeholder entityset
      EntitySet entitySet () const
      {
        assert(false && "Should not be used!");
        return EntitySet{};
      }

      /// Dummy function returning a placeholder entityset
      GlobalCoordinate operator() (GlobalCoordinate const&) const
      {
        assert(false && "Should not be used!");
        return GlobalCoordinate{};
      }

    private:
      /// All point coordinates inclusing the higher-order lagrange points
      Nodes nodes_;

      /// Parametrization for all elements
      Parametrization parametrization_;
    };

    // deduction guides
    template <class Grid>
    LagrangeGridCreator(GridFactory<Grid>&)
      -> LagrangeGridCreator<Grid>;

    template <class GridType, class FieldType, class Context>
    struct AssociatedGridFunction<LagrangeGridCreator<GridType>, FieldType, Context>
    {
      using type = LagrangeGridFunction<GridType, FieldType, Context>;
    };

    template <class Grid>
    class LagrangeGridCreator<Grid>::LocalParametrization
    {
      using ctype = typename Grid::ctype;

      using GlobalCoordinate = typename Grid::template Codim<0>::Entity::Geometry::GlobalCoordinate;
      using LocalCoordinate = typename Grid::template Codim<0>::Entity::Geometry::LocalCoordinate;
      using LocalGeometry = MultiLinearGeometry<ctype,Grid::dimension,Grid::dimension>;

      using LocalFE = LagrangeLocalFiniteElement<Vtk::LagrangePointSet, Grid::dimension, ctype, ctype>;
      using LocalBasis = typename LocalFE::Traits::LocalBasisType;
      using LocalBasisTraits = typename LocalBasis::Traits;

    public:
      /// Construct a local element parametrization
      template <class Nodes, class LocalParam>
      LocalParametrization (Nodes const& nodes, LocalParam const& param, int order)
        : localFE_(param.type, order)
        , localNodes_(param.nodes.size())
      {
        for (std::size_t i = 0; i < localNodes_.size(); ++i)
          localNodes_[i] = nodes[param.nodes[i]];
      }

      /// Construct a local element parametrization for elements with permuted corners
      template <class Nodes, class LocalParam, class LG>
      LocalParametrization (Nodes const& nodes, LocalParam const& param, int order, LG&& localGeometry)
        : LocalParametrization(nodes, param, order)
      {
        localGeometry_.emplace(std::forward<LG>(localGeometry));
      }

      /// Evaluate the local parametrization in local coordinates
      template <class LocalCoordinate>
      GlobalCoordinate operator() (LocalCoordinate const& local) const
      {
        // map coordinates if element corners are permuted
        LocalCoordinate x = localGeometry_ ? localGeometry_->global(local) : local;

        LocalBasis const& localBasis = localFE_.localBasis();
        localBasis.evaluateFunction(x, shapeValues_);
        assert(shapeValues_.size() == localNodes_.size());

        using field_type = typename LocalBasisTraits::RangeType::field_type;

        GlobalCoordinate out(0);
        for (std::size_t i = 0; i < shapeValues_.size(); ++i)
          out.axpy(field_type(shapeValues_[i]), localNodes_[i]);

        return out;
      }

    private:
      LocalFE localFE_;
      std::vector<GlobalCoordinate> localNodes_;
      std::optional<LocalGeometry> localGeometry_;

      mutable std::vector<typename LocalBasisTraits::RangeType> shapeValues_;
    };


    template <class Grid>
    class LagrangeGridCreator<Grid>::LocalFunction
    {
      using ctype = typename Grid::ctype;
      using LocalContext = typename Grid::template Codim<0>::Entity;
      using GlobalCoordinate = typename LocalContext::Geometry::GlobalCoordinate;
      using LocalCoordinate = typename LocalContext::Geometry::LocalCoordinate;
      using LocalParametrization = typename LagrangeGridCreator::LocalParametrization;

    public:
      explicit LocalFunction (LagrangeGridCreator const& gridCreator)
        : gridCreator_(&gridCreator)
      {}

      explicit LocalFunction (LagrangeGridCreator&& gridCreator) = delete;

      /// Collect a local parametrization on the element
      void bind (LocalContext const& element)
      {
        localContext_ = element;
        localParametrization_.emplace(gridCreator_->localParametrization(element));
      }

      void unbind () { /* do nothing */ }

      /// Evaluate the local parametrization in local coordinates
      GlobalCoordinate operator() (LocalCoordinate const& local) const
      {
        assert(!!localParametrization_);
        return (*localParametrization_)(local);
      }

      /// Return the bound element
      LocalContext const& localContext () const
      {
        return localContext_;
      }

    private:
      LagrangeGridCreator const* gridCreator_;

      LocalContext localContext_;
      std::optional<LocalParametrization> localParametrization_;
    };

  } // end namespace Vtk
} // end namespace Dune
