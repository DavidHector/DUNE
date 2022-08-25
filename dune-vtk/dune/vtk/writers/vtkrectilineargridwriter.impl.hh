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
void VtkRectilinearGridWriter<GV,DC>
  ::writeSerialFile (std::ofstream& out) const
{
  std::vector<pos_type> offsets; // pos => offset
  this->writeHeader(out, "RectilinearGrid");

  auto const& wholeExtent = dataCollector_->wholeExtent();
  out << "<RectilinearGrid"
      << " WholeExtent=\"" << Vtk::join(wholeExtent.begin(), wholeExtent.end()) << "\""
      << ">\n";

  dataCollector_->writeLocalPiece([&out](auto const& extent) {
    out << "<Piece Extent=\"" << Vtk::join(extent.begin(), extent.end()) << "\">\n";
  });

  // Write point coordinates for x, y, and z ordinate
  out << "<Coordinates>\n";
  writeCoordinates(out, offsets);
  out << "</Coordinates>\n";

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
  out << "</RectilinearGrid>\n";

  this->writeAppended(out, offsets);
  out << "</VTKFile>";
}


template <class GV, class DC>
void VtkRectilinearGridWriter<GV,DC>
  ::writeParallelFile (std::ofstream& out, std::string const& pfilename, int /*size*/) const
{
  this->writeHeader(out, "PRectilinearGrid");

  auto const& wholeExtent = dataCollector_->wholeExtent();
  out << "<PRectilinearGrid"
      << " GhostLevel=\"" << dataCollector_->ghostLevel() << "\""
      << " WholeExtent=\"" << Vtk::join(wholeExtent.begin(), wholeExtent.end()) << "\""
      << ">\n";

  // Write point coordinates for x, y, and z ordinate
  out << "<PCoordinates>\n";
  out << "<PDataArray Name=\"x\" type=\"" << to_string(datatype_) << "\" />\n";
  out << "<PDataArray Name=\"y\" type=\"" << to_string(datatype_) << "\" />\n";
  out << "<PDataArray Name=\"z\" type=\"" << to_string(datatype_) << "\" />\n";
  out << "</PCoordinates>\n";

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

  out << "</PRectilinearGrid>\n";
  out << "</VTKFile>";
}


template <class GV, class DC>
void VtkRectilinearGridWriter<GV,DC>
  ::writeCoordinates (std::ofstream& out, std::vector<pos_type>& offsets,
                      std::optional<std::size_t> timestep) const
{
  std::string names = "xyz";
  if (format_ == Vtk::FormatTypes::ASCII) {
    auto coordinates = dataCollector_->template coordinates<double>();
    for (std::size_t d = 0; d < 3; ++d) {
      out << "<DataArray type=\"" << to_string(datatype_) << "\" Name=\"" << names[d] << "\" format=\"ascii\"";
      if (timestep)
        out << " TimeStep=\"" << *timestep << "\"";
      out << ">\n";
      this->writeValuesAscii(out, coordinates[d]);
      out << "</DataArray>\n";
    }
  }
  else { // Vtk::FormatTypes::APPENDED format
    for (std::size_t j = 0; j < 3; ++j) {
      out << "<DataArray type=\"" << to_string(datatype_) << "\" Name=\"" << names[j] << "\" format=\"appended\"";
      if (timestep)
        out << " TimeStep=\"" << *timestep << "\"";
      out << " offset=";
      offsets.push_back(out.tellp());
      out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
      out << "/>\n";
    }
  }
}


template <class GV, class DC>
void VtkRectilinearGridWriter<GV,DC>
  ::writeGridAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const
{
  assert(is_a(format_, Vtk::FormatTypes::APPENDED) && "Function should by called only in appended mode!\n");

  // write coordinates along axis
  Vtk::mapDataTypes<std::is_floating_point, std::is_integral>(datatype_, headertype_,
  [&](auto f, auto h) {
    using F = typename decltype(f)::type;
    using H = typename decltype(h)::type;
    auto coordinates = dataCollector_->template coordinates<F>();
    blocks.push_back(this->template writeValuesAppended<H>(out, coordinates[0]));
    blocks.push_back(this->template writeValuesAppended<H>(out, coordinates[1]));
    blocks.push_back(this->template writeValuesAppended<H>(out, coordinates[2]));
  });
}

} // end namespace Dune
