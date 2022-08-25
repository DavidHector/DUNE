#pragma once

namespace Dune
{
  namespace Vtk
  {
    /// Context indicating that a GridFunction generates a local-function from point data
    struct PointContext {};

    /// Context indicating that a GridFunction generates a local-function from cell data
    struct CellContext {};

    /// \brief Type-Traits to associate a GridFunction to a GridCreator.
    /**
     * Each GridCreator type should specialize this template and set `type` to the
     * corresponding GridFunction type, e.g. Vtk::ContinuousGridFunction or
     * Vtk::LagrangeGridFunction.
     *
     * \tparam GridCreator   A Type implementing the GridCreatorInterface
     * \tparam FieldType     Coefficient type of the data extracted from the file.
     * \tparam Context       A context-type for specialization of the local-function, e.g.,
     *                       Vtk::PointContext or Vtk::CellContext.
     **/
    template <class GridCreator, class FieldType, class Context>
    struct AssociatedGridFunction
    {
      using type = void;
    };

  } // end namespace Vtk
} // end namespace Dune
