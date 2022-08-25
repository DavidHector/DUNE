#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/gridfactory.hh>

namespace Dune
{
  namespace Vtk
  {
    /// Base class for grid creators in a CRTP style.
    /**
    * Construct a grid from data read from VTK files.
    *
    * \tparam GridType         Model of Dune::Grid
    * \tparam GlobalCoordType  Type of the global coordinates.
    * \tparam DerivedType      Implementation of a concrete GridCreator.
    **/
    template <class GridType, class DerivedType>
    class GridCreatorInterface
    {
    public:
      using Grid = GridType;
      using GlobalCoordinate = typename Grid::template Codim<0>::Entity::Geometry::GlobalCoordinate;
      using Derived = DerivedType;

    public:
      /// Constructor. Stores a reference to the passed GridFactory.
      GridCreatorInterface (GridFactory<Grid>& factory)
        : factory_(stackobject_to_shared_ptr(factory))
      {}

      /// Constructor. Store the shared_ptr to the GridFactory.
      GridCreatorInterface (std::shared_ptr<GridFactory<Grid>> factory)
        : factory_(std::move(factory))
      {}

      /// Constructor. Construct a new GridFactory from the passed arguments.
      template <class... Args,
        std::enable_if_t<std::is_constructible<GridFactory<Grid>, Args...>::value,int> = 0>
      GridCreatorInterface (Args&&... args)
        : factory_(std::make_shared<GridFactory<Grid>>(std::forward<Args>(args)...))
      {}

      /// Insert all points as vertices into the factory
      void insertVertices (std::vector<GlobalCoordinate> const& points,
                          std::vector<std::uint64_t> const& point_ids)
      {
        asDerived().insertVerticesImpl(points, point_ids);
      }

      /// Create elements based on type and connectivity description
      void insertElements (std::vector<std::uint8_t> const& types,
                          std::vector<std::int64_t> const& offsets,
                          std::vector<std::int64_t> const& connectivity)
      {
        asDerived().insertElementsImpl(types, offsets, connectivity);
      }

      /// Insert part of a grid stored in file into factory
      void insertPieces (std::vector<std::string> const& pieces)
      {
        asDerived().insertPiecesImpl(pieces);
      }

      /// Construct the actual grid using the GridFactory
      std::unique_ptr<Grid> createGrid () const
      {
        return std::unique_ptr<Grid>(factory_->createGrid());
      }

      /// Return the associated GridFactory
      GridFactory<Grid>& factory ()
      {
        return *factory_;
      }

      /// Return the associated (const) GridFactory
      GridFactory<Grid> const& factory () const
      {
        return *factory_;
      }

      /// Return the mpi collective communicator
      auto comm () const
      {
        return MPIHelper::getCollectiveCommunication();
      }

    protected: // cast to derived type

      Derived& asDerived ()
      {
        return static_cast<Derived&>(*this);
      }

      const Derived& asDerived () const
      {
        return static_cast<const Derived&>(*this);
      }

    public: // default implementations

      void insertVerticesImpl (std::vector<GlobalCoordinate> const&,
                              std::vector<std::uint64_t> const&)
      {
        /* do nothing */
      }

      void insertElementsImpl (std::vector<std::uint8_t> const&,
                              std::vector<std::int64_t> const&,
                              std::vector<std::int64_t> const&)
      {
        /* do nothing */
      }

      void insertPiecesImpl (std::vector<std::string> const&)
      {
        /* do nothing */;
      }

    protected:
      std::shared_ptr<GridFactory<Grid>> factory_;
    };

  } // end namespace Vtk
} // end namespace Dune
