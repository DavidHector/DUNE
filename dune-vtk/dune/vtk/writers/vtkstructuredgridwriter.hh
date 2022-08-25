#pragma once

#include <array>
#include <iosfwd>
#include <map>

#include <dune/vtk/filewriter.hh>
#include <dune/vtk/function.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/datacollectors/structureddatacollector.hh>
#include <dune/vtk/utility/concepts.hh>

#include <dune/vtk/vtkwriterinterface.hh>

namespace Dune
{
  /// File-Writer for StructuredGrid VTK .vts files
  /**
   * Requirement:
   * - DataCollector must be a model of \ref StructuredDataCollector
   **/
  template <class GridView, class DataCollector = Vtk::StructuredDataCollector<GridView>>
  class VtkStructuredGridWriter
      : public VtkWriterInterface<GridView, DataCollector>
  {
    using Super = VtkWriterInterface<GridView, DataCollector>;
    using pos_type = typename Super::pos_type;

  public:
    /// forwarding constructor to \ref VtkWriterInterface
    using Super::Super;

  private:
    /// Write a serial VTK file in Unstructured format
    virtual void writeSerialFile (std::ofstream& out) const override;

    /// Write a parallel VTK file `pfilename.pvtu` in Unstructured format,
    /// with `size` the number of pieces and serial files given by `pfilename_p[i].vtu`
    /// for [i] in [0,...,size).
    virtual void writeParallelFile (std::ofstream& out, std::string const& pfilename, int size) const override;

    virtual std::string fileExtension () const override
    {
      return "vts";
    }

    virtual void writeGridAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const override;

  private:
    using Super::dataCollector_;
    using Super::format_;
    using Super::datatype_;
    using Super::headertype_;

    // attached data
    using Super::pointData_;
    using Super::cellData_;
  };

  // deduction guides
  template <class GridView, class... Args,
    Vtk::IsGridView<GridView> = true>
  VtkStructuredGridWriter(GridView, Args...)
    -> VtkStructuredGridWriter<GridView, Vtk::StructuredDataCollector<GridView>>;

  template <class DataCollector, class... Args,
    Vtk::IsDataCollector<DataCollector> = true>
  VtkStructuredGridWriter(DataCollector&, Args...)
    -> VtkStructuredGridWriter<typename DataCollector::GridView, DataCollector>;

  template <class DataCollector, class... Args,
    Vtk::IsDataCollector<DataCollector> = true>
  VtkStructuredGridWriter(std::shared_ptr<DataCollector>, Args...)
    -> VtkStructuredGridWriter<typename DataCollector::GridView, DataCollector>;

} // end namespace Dune

#include "vtkstructuredgridwriter.impl.hh"
