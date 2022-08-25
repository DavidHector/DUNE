#pragma once

#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/vtk/utility/enum.hh>
#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/string.hh>

namespace Dune {

template <class GV, class DC>
void VtkImageDataWriter<GV,DC>
  ::writeSerialFile (std::ofstream& out) const
{
  std::vector<pos_type> offsets; // pos => offset
  this->writeHeader(out, "ImageData");

  auto const& wholeExtent = dataCollector_->wholeExtent();
  auto const& origin = dataCollector_->origin();
  auto const& spacing = dataCollector_->spacing();
  out << "<ImageData"
      << " WholeExtent=\"" << Vtk::join(wholeExtent.begin(), wholeExtent.end()) << "\""
      << " Origin=\"" << Vtk::join(origin.begin(), origin.end()) << "\""
      << " Spacing=\"" << Vtk::join(spacing.begin(), spacing.end()) << "\""
      << ">\n";

  dataCollector_->writeLocalPiece([&out](auto const& extent) {
    out << "<Piece Extent=\"" << Vtk::join(extent.begin(), extent.end()) << "\">\n";
  });

  // Write data associated with grid points
  out << "<PointData" << this->getNames(pointData_) << ">\n";
  for (auto const& v : pointData_)
    this->writeData(out, offsets, v, Super::POINT_DATA);
  out << "</PointData>\n";

  // Write data associated with grid cells
  out << "<CellData" << this->getNames(cellData_) << ">\n";
  for (auto const& v : cellData_)
    this->writeData(out, offsets, v, Super::CELL_DATA);
  out << "</CellData>\n";

  out << "</Piece>\n";
  out << "</ImageData>\n";

  this->writeAppended(out, offsets);
  out << "</VTKFile>";
}


template <class GV, class DC>
void VtkImageDataWriter<GV,DC>
  ::writeParallelFile (std::ofstream& out, std::string const& pfilename, int /*size*/) const
{
  this->writeHeader(out, "PImageData");

  auto const& wholeExtent = dataCollector_->wholeExtent();
  auto const& origin = dataCollector_->origin();
  auto const& spacing = dataCollector_->spacing();
  out << "<PImageData"
      << " GhostLevel=\"" << dataCollector_->ghostLevel() << "\""
      << " WholeExtent=\"" << Vtk::join(wholeExtent.begin(), wholeExtent.end()) << "\""
      << " Origin=\"" << Vtk::join(origin.begin(), origin.end()) << "\""
      << " Spacing=\"" << Vtk::join(spacing.begin(), spacing.end()) << "\""
      << ">\n";

  // Write data associated with grid points
  out << "<PPointData" << this->getNames(pointData_) << ">\n";
  for (auto const& v : pointData_) {
    out << "<PDataArray"
        << " Name=\"" << v.name() << "\""
        << " type=\"" << to_string(v.dataType()) << "\""
        << " NumberOfComponents=\"" << v.numComponents() << "\""
        << " />\n";
  }
  out << "</PPointData>\n";

  // Write data associated with grid cells
  out << "<PCellData" << this->getNames(cellData_) << ">\n";
  for (auto const& v : cellData_) {
    out << "<PDataArray"
        << " Name=\"" << v.name() << "\""
        << " type=\"" <<  to_string(v.dataType()) << "\""
        << " NumberOfComponents=\"" << v.numComponents() << "\""
        << " />\n";
  }
  out << "</PCellData>\n";

  // Write piece file references
  dataCollector_->writePieces([&out,pfilename,ext=this->fileExtension()](int p, auto const& extent, bool write_extent)
  {
    std::string piece_source = pfilename + "_p" + std::to_string(p) + "." + ext;
    out << "<Piece Source=\"" << piece_source << "\"";
    if (write_extent)
      out << " Extent=\"" << Vtk::join(extent.begin(), extent.end()) << "\"";
     out << " />\n";
  });

  out << "</PImageData>\n";
  out << "</VTKFile>";
}

} // end namespace Dune
