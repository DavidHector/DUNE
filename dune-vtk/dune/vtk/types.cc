#include <config.h>

#include <iostream>
#include <map>

#include <dune/common/exceptions.hh>
#include <dune/vtk/types.hh>

namespace Dune {
namespace Vtk {

std::string to_string (Vtk::FormatTypes type)
{
  switch (type) {
    case Vtk::FormatTypes::ASCII:      return "ascii";
    case Vtk::FormatTypes::BINARY:     return "binary";
    case Vtk::FormatTypes::COMPRESSED: return "compressed";
    case Vtk::FormatTypes::APPENDED:   return "appended";
    default:
      DUNE_THROW(RangeError, "FormatType not found.");
      std::abort();
  }
}

Vtk::FormatTypes formatTypeOf (Dune::VTK::OutputType o)
{
  switch (o) {
    case Dune::VTK::ascii:          return Vtk::FormatTypes::ASCII;
 // case Dune::VTK::base64:         return Vtk::FormatTypes::BASE64;
    case Dune::VTK::appendedraw:    return Vtk::FormatTypes::BINARY;
 // case Dune::VTK::appendedbase64: return Vtk::FormatTypes::BASE64;
 // case Dune::VTK::binarycompressed:   return Vtk::FormatTypes::COMPRESSED;
 // case Dune::VTK::compressedappended: return Vtk::FormatTypes::COMPRESSED;
    default:
      DUNE_THROW(RangeError, "OutputType not supported.");
      std::abort();
  }
}


std::string to_string (Vtk::RangeTypes type)
{
  switch (type) {
    case Vtk::RangeTypes::UNSPECIFIED:  return "unspecified";
    case Vtk::RangeTypes::AUTO:         return "auto";
    case Vtk::RangeTypes::SCALAR:       return "scalar";
    case Vtk::RangeTypes::VECTOR:       return "vector";
    case Vtk::RangeTypes::TENSOR:       return "tensor";
    default:
      DUNE_THROW(RangeError, "RangeType not found.");
      std::abort();
  }
}

// Map a dune-grid FieldInfo::Type to ValueTypes
Vtk::RangeTypes rangeTypeOf (Dune::VTK::FieldInfo::Type t)
{
  switch (t) {
    case Dune::VTK::FieldInfo::Type::scalar: return Vtk::RangeTypes::UNSPECIFIED;
    case Dune::VTK::FieldInfo::Type::vector: return Vtk::RangeTypes::VECTOR;
    case Dune::VTK::FieldInfo::Type::tensor: return Vtk::RangeTypes::TENSOR;
    default:
      DUNE_THROW(RangeError, "FieldInfo::Type not supported.");
      std::abort();
  }
}

// Map a number of components to a corresponding value type
Vtk::RangeTypes rangeTypeOf (int ncomps)
{
  return ncomps > 9 ? Vtk::RangeTypes::UNSPECIFIED :
         ncomps > 3 ? Vtk::RangeTypes::TENSOR :
         ncomps > 1 ? Vtk::RangeTypes::VECTOR :
                      Vtk::RangeTypes::SCALAR;
}


std::string to_string (Vtk::DataTypes type)
{
  switch (type) {
    case Vtk::DataTypes::UNKNOWN: return "unknown";
    case Vtk::DataTypes::INT8:    return "Int8";
    case Vtk::DataTypes::UINT8:   return "UInt8";
    case Vtk::DataTypes::INT16:   return "Int16";
    case Vtk::DataTypes::UINT16:  return "UInt16";
    case Vtk::DataTypes::INT32:   return "Int32";
    case Vtk::DataTypes::UINT32:  return "UInt32";
    case Vtk::DataTypes::INT64:   return "Int64";
    case Vtk::DataTypes::UINT64:  return "UInt64";
    case Vtk::DataTypes::FLOAT32: return "Float32";
    case Vtk::DataTypes::FLOAT64: return "Float64";
    default:
      DUNE_THROW(RangeError, "DataType not found.");
      std::abort();
  }
}


Vtk::DataTypes dataTypeOf (Dune::VTK::Precision p)
{
  switch (p) {
    case Dune::VTK::Precision::int32:    return Vtk::DataTypes::INT32;
    case Dune::VTK::Precision::uint8:    return Vtk::DataTypes::UINT8;
    case Dune::VTK::Precision::uint32:   return Vtk::DataTypes::UINT32;
    case Dune::VTK::Precision::float32:  return Vtk::DataTypes::FLOAT32;
    case Dune::VTK::Precision::float64:  return Vtk::DataTypes::FLOAT64;
    default:
      DUNE_THROW(RangeError, "Precision not supported.");
      std::abort();
  }
}


Vtk::DataTypes dataTypeOf (std::string s)
{
  static const std::map<std::string, Vtk::DataTypes> to_datatype{
    {"Int8",    Vtk::DataTypes::INT8},
    {"UInt8",   Vtk::DataTypes::UINT8},
    {"Int16",   Vtk::DataTypes::INT16},
    {"UInt16",  Vtk::DataTypes::UINT16},
    {"Int32",   Vtk::DataTypes::INT32},
    {"UInt32",  Vtk::DataTypes::UINT32},
    {"Int64",   Vtk::DataTypes::INT64},
    {"UInt64",  Vtk::DataTypes::UINT64},
    {"Float32", Vtk::DataTypes::FLOAT32},
    {"Float64", Vtk::DataTypes::FLOAT64}
  };
  auto it = to_datatype.find(s);
  return it != to_datatype.end() ? it->second : Vtk::DataTypes::UNKNOWN;
}


std::string to_string (Vtk::CompressorTypes type)
{
  switch (type) {
    case Vtk::CompressorTypes::ZLIB: return "vtkZLibDataCompressor";
    case Vtk::CompressorTypes::LZ4:  return "vtkLZ4DataCompressor";
    case Vtk::CompressorTypes::LZMA: return "vtkLZMADataCompressor";
    default:
      DUNE_THROW(RangeError, "CompressorTypes not found.");
      std::abort();
  }
}


GeometryType to_geometry (std::uint8_t cell)
{
  switch (cell) {
    case Vtk::CellType::VERTEX:     return GeometryTypes::vertex;
    case Vtk::CellType::LINE:       return GeometryTypes::line;
    case Vtk::CellType::TRIANGLE:   return GeometryTypes::triangle;
    case Vtk::CellType::QUAD:       return GeometryTypes::quadrilateral;
    case Vtk::CellType::TETRA:      return GeometryTypes::tetrahedron;
    case Vtk::CellType::HEXAHEDRON: return GeometryTypes::hexahedron;
    case Vtk::CellType::WEDGE:      return GeometryTypes::prism;
    case Vtk::CellType::PYRAMID:    return GeometryTypes::pyramid;

    // Quadratic VTK cell types
    case Vtk::CellType::QUADRATIC_EDGE:        return GeometryTypes::line;
    case Vtk::CellType::QUADRATIC_TRIANGLE:    return GeometryTypes::triangle;
    case Vtk::CellType::QUADRATIC_QUAD:        return GeometryTypes::quadrilateral;
    case Vtk::CellType::QUADRATIC_TETRA:       return GeometryTypes::tetrahedron;
    case Vtk::CellType::QUADRATIC_HEXAHEDRON:  return GeometryTypes::hexahedron;

    // Arbitrary order Lagrange elements
    case Vtk::CellType::LAGRANGE_CURVE:        return GeometryTypes::line;
    case Vtk::CellType::LAGRANGE_TRIANGLE:     return GeometryTypes::triangle;
    case Vtk::CellType::LAGRANGE_QUADRILATERAL:return GeometryTypes::quadrilateral;
    case Vtk::CellType::LAGRANGE_TETRAHEDRON:  return GeometryTypes::tetrahedron;
    case Vtk::CellType::LAGRANGE_HEXAHEDRON:   return GeometryTypes::hexahedron;
    case Vtk::CellType::LAGRANGE_WEDGE:        return GeometryTypes::prism;

    default:
      DUNE_THROW(RangeError, "CellType does not map to GeometryType.");
      std::abort();
  }
}


CellType::CellType (GeometryType const& t, Parametrization parametrization)
  : noPermutation_(true)
{
  if (parametrization == LINEAR) {
    if (t.isVertex()) {
      type_ = VERTEX;
      permutation_ = {0};
    }
    else if (t.isLine()) {
      type_ = LINE;
      permutation_ = {0,1};
    }
    else if (t.isTriangle()) {
      type_ = TRIANGLE;
      permutation_ = {0,1,2};
    }
    else if (t.isQuadrilateral()) {
      type_ = QUAD;
      permutation_ = {0,1,3,2};
      noPermutation_ = false;
    }
    else if (t.isTetrahedron()) {
      type_ = TETRA;
      permutation_ = {0,1,2,3};
    }
    else if (t.isHexahedron()) {
      type_ = HEXAHEDRON;
      permutation_ = {0,1,3,2,4,5,7,6};
      noPermutation_ = false;
    }
    else if (t.isPrism()) {
      type_ = WEDGE;
      permutation_ = {0,2,1,3,5,4};
      noPermutation_ = false;
    }
    else if (t.isPyramid()) {
      type_ = PYRAMID;
      permutation_ = {0,1,3,2,4};
      noPermutation_ = false;
    }
    else if (t.isNone() && t.dim() == 1) {
      type_ = LINE;
      permutation_ = {0,1};
    }
    else if (t.isNone() && t.dim() == 2) {
      type_ = POLYGON;
      permutation_ = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    }
    else {
      std::cerr << "Geometry Type not supported by VTK!\n";
      std::abort();
    }
  } else if (parametrization == QUADRATIC) {
    if (t.isLine()) {
      type_ = QUADRATIC_EDGE;
      permutation_ = {0,1, 0};
    }
    else if (t.isTriangle()) {
      type_ = QUADRATIC_TRIANGLE;
      permutation_ = {0,1,2, 0,2,1};
      noPermutation_ = false;
    }
    else if (t.isQuadrilateral()) {
      type_ = QUADRATIC_QUAD;
      permutation_ = {0,1,3,2, 2,1,3,0};
      noPermutation_ = false;
    }
    else if (t.isTetrahedron()) {
      type_ = QUADRATIC_TETRA;
      permutation_ = {0,1,2,3, 0,2,1,3,4,5};
      noPermutation_ = false;
    }
    else if (t.isHexahedron()) {
      type_ = QUADRATIC_HEXAHEDRON;
      permutation_ = {0,1,3,2,4,5,7,6, 6,5,7,4,10,9,11,8,0,1,3,2};
      noPermutation_ = false;
    }
    else {
      std::cerr << "Geometry Type not supported by VTK!\n";
      std::abort();
    }
  } else if (parametrization == LAGRANGE) {
    if (t.isLine()) {
      type_ = LAGRANGE_CURVE;
    }
    else if (t.isTriangle()) {
      type_ = LAGRANGE_TRIANGLE;
    }
    else if (t.isQuadrilateral()) {
      type_ = LAGRANGE_QUADRILATERAL;
    }
    else if (t.isTetrahedron()) {
      type_ = LAGRANGE_TETRAHEDRON;
    }
    else if (t.isHexahedron()) {
      type_ = LAGRANGE_HEXAHEDRON;
    }
    else if (t.isPrism()) {
      type_ = LAGRANGE_WEDGE;
    }
    else if (t.isPyramid()) {
      type_ = LAGRANGE_PYRAMID;
    }
    else {
      std::cerr << "Geometry Type not supported by VTK!\n";
      std::abort();
    }
  }
}

} } // end namespace Dune::Vtk
