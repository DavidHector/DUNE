#pragma once

#include <cassert>
#include <array>

#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>

namespace Dune
{
  namespace Vtk
  {
    namespace Impl
    {
      // forward declaration
      template <class K, unsigned int dim>
      class LagrangePointSetBuilder;
    }


    /// \brief A set of lagrange points compatible with the numbering of VTK and Gmsh
    /**
     * \tparam K    Field-type for the coordinates
     * \tparam dim  Dimension of the coordinates
     **/
    template <class K, unsigned int dim>
    class LagrangePointSet
        : public EmptyPointSet<K, dim>
    {
      using Super = EmptyPointSet<K, dim>;

    public:
      static const unsigned int dimension = dim;

      LagrangePointSet (std::size_t order)
        : Super(order)
      {
        assert(order > 0);
      }

      /// Fill the lagrange points for the given geometry type
      void build (GeometryType gt)
      {
        assert(gt.dim() == dimension);
        builder_(gt, order(), points_);
      }

      /// Fill the lagrange points for the given topology type `Topology`
      template <GeometryType::Id geometryId>
      bool build ()
      {
        build(GeometryType(geometryId));
        return true;
      }

      /// Returns whether the point set support the given topology type `Topology` and can
      /// generate point for the given order.
      template <GeometryType::Id geometryId>
      static bool supports (std::size_t order)
      {
        return true;
      }

      using Super::order;

    private:
      using Super::points_;
      Impl::LagrangePointSetBuilder<K,dim> builder_;
    };


    namespace Impl
    {
      // Build for lagrange point sets in different dimensions
      // Specialized for dim=1,2,3
      template <class K, unsigned int dim>
      class LagrangePointSetBuilder
      {
      public:
        template <class Points>
        void operator()(GeometryType, unsigned int, Points& points) const
        {
          DUNE_THROW(Dune::NotImplemented,
            "Lagrange points not yet implemented for this GeometryType.");
        }
      };


      // Lagrange points on point geometries
      template <class K>
      class LagrangePointSetBuilder<K,0>
      {
        static constexpr int dim = 0;
        using LP = LagrangePoint<K,dim>;
        using Vec = typename LP::Vector;
        using Key = LocalKey;

      public:
        template <class Points>
        void operator()(GeometryType gt, int /*order*/, Points& points) const;
      };


      // Lagrange points on line geometries
      template <class K>
      class LagrangePointSetBuilder<K,1>
      {
        static constexpr int dim = 1;
        using LP = LagrangePoint<K,dim>;
        using Vec = typename LP::Vector;
        using Key = LocalKey;

      public:
        template <class Points>
        void operator()(GeometryType gt, int order, Points& points) const;
      };


      // Lagrange points on 2d geometries
      template <class K>
      class LagrangePointSetBuilder<K,2>
      {
        static constexpr int dim = 2;
        using LP = LagrangePoint<K,dim>;
        using Vec = typename LP::Vector;
        using Key = LocalKey;

        friend class LagrangePointSetBuilder<K,3>;

      public:
        template <class Points>
        void operator()(GeometryType gt, int order, Points& points) const;

      private: // implementation details

        // Construct the point set in a triangle element.
        // Loop from the outside to the inside
        template <class Points>
        void buildTriangle (std::size_t nPoints, int order, Points& points) const;

        // "Barycentric index" is a triplet of integers, each running from 0 to
        // <Order>. It is the index of a point on the triangle in barycentric
        // coordinates.
        static void barycentricIndex (int index, std::array<int,3>& bindex, int order);

        // Construct the point set in the quad element
        // 1. build equispaced points with index tuple (i,j)
        // 2. map index tuple to DOF index and LocalKey
        template <class Points>
        void buildQuad(std::size_t nPoints, int order, Points& points) const;

        // Obtain the VTK DOF index of the node (i,j) in the quad element
        // and construct a LocalKey
        static std::pair<int,Key> calcQuadKey (int i, int j, std::array<int,2> order);
      };


      // Lagrange points on 3d geometries
      template <class K>
      class LagrangePointSetBuilder<K,3>
      {
        static constexpr int dim = 3;
        using LP = LagrangePoint<K,dim>;
        using Vec = typename LP::Vector;
        using Key = LocalKey;

      public:
        template <class Points>
        void operator() (GeometryType gt, unsigned int order, Points& points) const;

      private: // implementation details

        // Construct the point set in the tetrahedron element
        // 1. construct barycentric (index) coordinates
        // 2. obtains the DOF index, LocalKey and actual coordinate from barycentric index
        template <class Points>
        void buildTetra (std::size_t nPoints, int order, Points& points) const;

        // "Barycentric index" is a set of 4 integers, each running from 0 to
        // <Order>. It is the index of a point in the tetrahedron in barycentric
        // coordinates.
        static void barycentricIndex (int p, std::array<int,4>& bindex, int order);

        // Construct the point set in the heyhedral element
        // 1. build equispaced points with index tuple (i,j,k)
        // 2. map index tuple to DOF index and LocalKey
        template <class Points>
        void buildHex (std::size_t nPoints, int order, Points& points) const;

        // Obtain the VTK DOF index of the node (i,j,k) in the hexahedral element
        static std::pair<int,Key> calcHexKey (int i, int j, int k, std::array<int,3> order);
      };

    } // end namespace Impl
  } // end namespace Vtk
} // end namespace Dune

#include <dune/vtk/utility/lagrangepoints.impl.hh>
