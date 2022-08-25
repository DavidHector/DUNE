#include <fstream>
#include <iterator>
#include <sstream>
#include <string>

#if HAVE_VTK_ZLIB
#include <zlib.h>
#endif

#include <dune/common/classname.hh>
#include <dune/common/version.hh>

#include "utility/errors.hh"
#include "utility/filesystem.hh"
#include "utility/string.hh"

namespace Dune {

template <class Grid, class Creator, class Field>
void VtkReader<Grid,Creator,Field>::read (std::string const& filename, bool fillCreator)
{
  // check whether file exists!
  if (!Vtk::exists(filename))
    DUNE_THROW(IOError, "File " << filename << " does not exist!");

  std::ifstream input(filename, std::ios_base::in | std::ios_base::binary);
  VTK_ASSERT(input.is_open());

  std::string ext = Vtk::Path(filename).extension().string();
  if (ext == ".vtu") {
    readSerialFileFromStream(input, fillCreator);
    pieces_.push_back(filename);
  } else if (ext == ".pvtu") {
    readParallelFileFromStream(input, comm().rank(), comm().size(), fillCreator);
  } else {
    DUNE_THROW(Dune::VtkError, "File has unknown file-extension '" << ext << "'. Allowed are only '.vtu' and '.pvtu'.");
  }
  read_ = true;
}


template <class Grid, class Creator, class Field>
void VtkReader<Grid,Creator,Field>::readSerialFileFromStream (std::ifstream& input, bool fillCreator)
{
  clear();
  compressor_ = Vtk::CompressorTypes::NONE;
  Vtk::DataTypes header_type = Vtk::DataTypes::UINT32;
  std::string data_id = "";
  std::string data_format = "";
  Vtk::DataTypes data_type = Vtk::DataTypes::UNKNOWN;
  unsigned int data_components = 0;
  std::uint64_t data_offset = 0;

  Sections section = NO_SECTION;
  for (std::string line; std::getline(input, line); ) {
    Vtk::ltrim(line);

    if (isSection(line, "VTKFile", section)) {
      bool closed = false;
      auto attr = parseXml(line, closed);

      if (!attr["type"].empty())
        VTK_ASSERT_MSG(attr["type"] == "UnstructuredGrid", "VtkReader supports UnstructuredGrid types");
      if (!attr["byte_order"].empty())
        VTK_ASSERT_MSG(attr["byte_order"] == "LittleEndian", "LittleEndian byte order supported");

      if (attr["header_type"] == "UInt32")
        header_type = Vtk::DataTypes::UINT32;
      else if (attr["header_type"] == "UInt64")
        header_type = Vtk::DataTypes::UINT64;

      if (attr["compressor"] == "vtkZLibDataCompressor")
        compressor_ = Vtk::CompressorTypes::ZLIB;
      else if (attr["compressor"] == "vtkLZ4DataCompressor")
        compressor_ = Vtk::CompressorTypes::LZ4;
      else if (attr["compressor"] == "vtkLZMADataCompressor")
        compressor_ = Vtk::CompressorTypes::LZMA;

      section = VTK_FILE;
    }
    else if (isSection(line, "/VTKFile", section, VTK_FILE)) {
      section = NO_SECTION;
      break;
    }
    else if (isSection(line, "UnstructuredGrid", section, VTK_FILE))
      section = UNSTRUCTURED_GRID;
    else if (isSection(line, "/UnstructuredGrid", section, UNSTRUCTURED_GRID))
      section = VTK_FILE;
    else if (isSection(line, "Piece", section, UNSTRUCTURED_GRID)) {
      bool closed = false;
      auto attr = parseXml(line, closed);

      VTK_ASSERT_MSG(attr.count("NumberOfPoints") > 0 && attr.count("NumberOfCells") > 0,
        "Number of points or cells in file must be > 0");
      numberOfPoints_ = std::stoul(attr["NumberOfPoints"]);
      numberOfCells_ = std::stoul(attr["NumberOfCells"]);
      section = PIECE;
    }
    else if (isSection(line, "/Piece", section, PIECE))
      section = UNSTRUCTURED_GRID;
    else if (isSection(line, "PointData", section, PIECE))
      section = POINT_DATA;
    else if (isSection(line, "/PointData", section, POINT_DATA))
      section = PIECE;
    else if (isSection(line, "CellData", section, PIECE))
      section = CELL_DATA;
    else if (isSection(line, "/CellData", section, CELL_DATA))
      section = PIECE;
    else if (isSection(line, "Points", section, PIECE))
      section = POINTS;
    else if (isSection(line, "/Points", section, POINTS))
      section = PIECE;
    else if (isSection(line, "Cells", section, PIECE))
      section = CELLS;
    else if (isSection(line, "/Cells", section, CELLS))
      section = PIECE;
    else if (line.substr(1,9) == "DataArray") {
      bool closed = false;
      auto attr = parseXml(line, closed);

      data_type = Vtk::dataTypeOf(attr["type"]);

      // Use Section.Name as id
      data_id = toString(section) + "." + attr["Name"];

      if (section == POINTS)
        // In the Points section must only be one DataArray with id=Points
        data_id = "Points";

      data_components = 1;
      if (!attr["NumberOfComponents"].empty())
        data_components = std::stoul(attr["NumberOfComponents"]);

      // determine FormatType
      data_format = Vtk::to_lower(attr["format"]);
      if (data_format == "appended") {
        format_ = compressor_ != Vtk::CompressorTypes::NONE ? Vtk::FormatTypes::COMPRESSED : Vtk::FormatTypes::BINARY;
      } else {
        format_ = Vtk::FormatTypes::ASCII;
      }

      // Offset makes sense in appended mode only
      data_offset = 0;
      if (!attr["offset"].empty()) {
        data_offset = std::stoul(attr["offset"]);
        VTK_ASSERT_MSG(data_format == "appended", "Attribute 'offset' only supported by appended mode");
      }

      // Store attributes of DataArray
      dataArray_[data_id] = {attr["Name"], data_type, data_components, data_offset, section};

      // Skip section in appended mode
      if (data_format == "appended") {
        if (!closed) {
          while (std::getline(input, line)) {
            Vtk::ltrim(line);
            if (line.substr(1,10) == "/DataArray")
              break;
          }
        }
        continue;
      }

      if (section == POINT_DATA)
        section = PD_DATA_ARRAY;
      else if (section == POINTS)
        section = POINTS_DATA_ARRAY;
      else if (section == CELL_DATA)
        section = CD_DATA_ARRAY;
      else if (section == CELLS)
        section = CELLS_DATA_ARRAY;
      else
        DUNE_THROW(Dune::VtkError, "Wrong section for <DataArray>");
    }
    else if (line.substr(1,10) == "/DataArray") {
      if (section == PD_DATA_ARRAY)
        section = POINT_DATA;
      else if (section == POINTS_DATA_ARRAY)
        section = POINTS;
      else if (section == CD_DATA_ARRAY)
        section = CELL_DATA;
      else if (section == CELLS_DATA_ARRAY)
        section = CELLS;
      else
        DUNE_THROW(Dune::VtkError, "Wrong section for </DataArray>");
    }
    else if (isSection(line, "AppendedData", section, VTK_FILE)) {
      bool closed = false;
      auto attr = parseXml(line, closed);
      if (!attr["encoding"].empty())
        VTK_ASSERT_MSG(attr["encoding"] == "raw", "base64 encoding not supported");

      offset0_ = findAppendedDataPosition(input);
      Vtk::mapDataTypes<std::is_floating_point, std::is_integral>(dataArray_["Points"].type, header_type,
      [&](auto f, auto h) {
        this->readPointsAppended(f,h,input);
        this->readCellsAppended(h,input);
      });

      // read point and cell data
      for (auto const& d : dataArray_) {
        if (d.second.section == POINT_DATA) {
          Vtk::mapDataTypes<std::is_floating_point, std::is_integral>(d.second.type, header_type,
          [&](auto f, auto h) {
            this->readPointDataAppended(f,h,input,d.first);
          });
        }
        else if (d.second.section == CELL_DATA) {
          Vtk::mapDataTypes<std::is_floating_point, std::is_integral>(d.second.type, header_type,
          [&](auto f, auto h) {
            this->readCellDataAppended(f,h,input,d.first);
          });
        }
      }

      section = NO_SECTION; // finish reading after appended section
      break;
    }
    else if (isSection(line, "/AppendedData", section, APPENDED_DATA))
      section = VTK_FILE;

    switch (section) {
      case PD_DATA_ARRAY:
        section = readPointData(input, data_id);
        break;
      case POINTS_DATA_ARRAY:
        section = readPoints(input, data_id);
        break;
      case CD_DATA_ARRAY:
        section = readCellData(input, data_id);
        break;
      case CELLS_DATA_ARRAY:
        section = readCells(input, data_id);
        break;
      default:
        break;
    }
  }

  if (section != NO_SECTION)
    DUNE_THROW(Dune::VtkError, "VTK-File is incomplete. It must end with </VTKFile>!");

  if (fillCreator)
    fillGridCreator();
}


template <class Grid, class Creator, class Field>
void VtkReader<Grid,Creator,Field>::readParallelFileFromStream (std::ifstream& input, int /* commRank */, int /* commSize */, bool fillCreator)
{
  clear();

  [[maybe_unused]] Vtk::DataTypes header_type = Vtk::DataTypes::UINT32;
  compressor_ = Vtk::CompressorTypes::NONE;

  Sections section = NO_SECTION;
  for (std::string line; std::getline(input, line); ) {
    Vtk::ltrim(line);

    if (isSection(line, "VTKFile", section)) {
      bool closed = false;
      auto attr = parseXml(line, closed);

      if (!attr["type"].empty())
        VTK_ASSERT_MSG(attr["type"] == "PUnstructuredGrid", "VtkReader supports PUnstructuredGrid types");
      if (!attr["version"].empty())
        VTK_ASSERT_MSG(std::stod(attr["version"]) == 1.0, "File format must be 1.0");
      if (!attr["byte_order"].empty())
        VTK_ASSERT_MSG(attr["byte_order"] == "LittleEndian", "LittleEndian byte order supported");

      if (attr["header_type"] == "UInt32")
        header_type = Vtk::DataTypes::UINT32;
      else if (attr["header_type"] == "UInt64")
        header_type = Vtk::DataTypes::UINT64;

      if (attr["compressor"] == "vtkZLibDataCompressor")
        compressor_ = Vtk::CompressorTypes::ZLIB;
      else if (attr["compressor"] == "vtkLZ4DataCompressor")
        compressor_ = Vtk::CompressorTypes::LZ4;
      else if (attr["compressor"] == "vtkLZMADataCompressor")
        compressor_ = Vtk::CompressorTypes::LZMA;

      section = VTK_FILE;
    }
    else if (isSection(line, "/VTKFile", section, VTK_FILE)) {
      section = NO_SECTION;
      break;
    } else if (isSection(line, "PUnstructuredGrid", section, VTK_FILE))
      section = UNSTRUCTURED_GRID;
    else if (isSection(line, "/PUnstructuredGrid", section, UNSTRUCTURED_GRID))
      section = VTK_FILE;
    else if (isSection(line, "Piece", section, UNSTRUCTURED_GRID)) {
      bool closed = false;
      auto attr = parseXml(line, closed);

      VTK_ASSERT_MSG(attr.count("Source") > 0, "No source files for partitions provided");
      pieces_.push_back(attr["Source"]);
    }
  }

  VTK_ASSERT_MSG(section == NO_SECTION, "VTK-File is incomplete. It must end with </VTKFile>!");

  if (fillCreator)
    fillGridCreator();
}


// @{ implementation detail
/**
 * Read ASCII data from `input` stream into vector `values`
 * \param max_size  Upper bound for the number of values
 * \param section   Current XML section you are reading in
 * \param parent_section   XML Section to return when current `section` is finished.
 **/
template <class IStream, class T, class Sections>
Sections readDataArray (IStream& input, std::vector<T>& values, std::size_t max_size,
                        Sections section, Sections parent_section)
{
  values.reserve(max_size < std::size_t(-1) ? max_size : 0);
  using S = std::conditional_t<(sizeof(T) <= 1), std::uint16_t, T>; // problem when reading chars as ints

  std::size_t idx = 0;
  for (std::string line; std::getline(input, line);) {
    Vtk::trim(line);
    if (line.substr(1,10) == "/DataArray")
      return parent_section;
    if (line[0] == '<')
      break;

    std::istringstream stream(line);
    S value;
    for (; stream >> value; idx++)
      values.push_back(T(value));
    if (idx >= max_size)
      break;
  }

  return section;
}

template <class IStream, class Sections>
Sections skipRestOfDataArray (IStream& input, Sections section, Sections parent_section)
{
  for (std::string line; std::getline(input, line);) {
    Vtk::ltrim(line);
    if (line.substr(1,10) == "/DataArray")
      return parent_section;
  }

  return section;
}
// @}


// Read values stored on the cells with name `name`
template <class Grid, class Creator, class Field>
typename VtkReader<Grid,Creator,Field>::Sections
VtkReader<Grid,Creator,Field>::readCellData (std::ifstream& input, std::string id)
{
  VTK_ASSERT(numberOfCells_ > 0);
  unsigned int components = dataArray_[id].components;

  Sections sec;
  std::vector<Field>& values = cellData_[id];
  sec = readDataArray(input, values, components*numberOfCells_, CD_DATA_ARRAY, CELL_DATA);
  if (sec != CELL_DATA)
    sec = skipRestOfDataArray(input, CD_DATA_ARRAY, CELL_DATA);
  VTK_ASSERT(sec == CELL_DATA);
  VTK_ASSERT(values.size() == components*numberOfCells_);

  return sec;
}


template <class Grid, class Creator, class Field>
  template <class FloatType, class HeaderType>
void VtkReader<Grid,Creator,Field>::readCellDataAppended (MetaType<FloatType>, MetaType<HeaderType>, std::ifstream& input, std::string id)
{
  VTK_ASSERT(numberOfCells_ > 0);
  unsigned int components = dataArray_[id].components;

  std::vector<FloatType> values;
  readAppended(input, values, HeaderType(dataArray_[id].offset));
  VTK_ASSERT(values.size() == components*numberOfCells_);

  cellData_[id].resize(values.size());
  std::copy(values.begin(), values.end(), cellData_[id].begin());
}


template <class Grid, class Creator, class Field>
typename VtkReader<Grid,Creator,Field>::Sections
VtkReader<Grid,Creator,Field>::readPointData (std::ifstream& input, std::string id)
{
  VTK_ASSERT(numberOfPoints_ > 0);
  unsigned int components = dataArray_[id].components;

  Sections sec;
  std::vector<Field>& values = pointData_[id];
  sec = readDataArray(input, values, components*numberOfPoints_, PD_DATA_ARRAY, POINT_DATA);
  if (sec != POINT_DATA)
    sec = skipRestOfDataArray(input, PD_DATA_ARRAY, POINT_DATA);
  VTK_ASSERT(sec == POINT_DATA);
  VTK_ASSERT(values.size() == components*numberOfPoints_);

  return sec;
}


template <class Grid, class Creator, class Field>
  template <class FloatType, class HeaderType>
void VtkReader<Grid,Creator,Field>::readPointDataAppended (MetaType<FloatType>, MetaType<HeaderType>, std::ifstream& input, std::string id)
{
  VTK_ASSERT(numberOfPoints_ > 0);
  unsigned int components = dataArray_[id].components;

  std::vector<FloatType> values;
  readAppended(input, values, HeaderType(dataArray_[id].offset));
  VTK_ASSERT(values.size() == components*numberOfPoints_);

  pointData_[id].resize(values.size());
  std::copy(values.begin(), values.end(), pointData_[id].begin());
}


template <class Grid, class Creator, class Field>
typename VtkReader<Grid,Creator,Field>::Sections
VtkReader<Grid,Creator,Field>::readPoints (std::ifstream& input, std::string id)
{
  using T = typename GlobalCoordinate::value_type;
  VTK_ASSERT(numberOfPoints_ > 0);
  VTK_ASSERT(id == "Points");
  VTK_ASSERT(dataArray_["Points"].components == 3u);

  Sections sec;

  std::vector<T> point_values;
  sec = readDataArray(input, point_values, 3*numberOfPoints_, POINTS_DATA_ARRAY, POINTS);
  if (sec != POINTS)
    sec = skipRestOfDataArray(input, POINTS_DATA_ARRAY, POINTS);
  VTK_ASSERT(sec == POINTS);
  VTK_ASSERT(point_values.size() == 3*numberOfPoints_);

  // extract points from continuous values
  GlobalCoordinate p;
  vec_points.reserve(numberOfPoints_);
  std::size_t idx = 0;
  for (std::size_t i = 0; i < numberOfPoints_; ++i) {
    for (std::size_t j = 0; j < p.size(); ++j)
      p[j] = point_values[idx++];
    idx += (3u - p.size());
    vec_points.push_back(p);
  }

  return sec;
}


template <class Grid, class Creator, class Field>
  template <class FloatType, class HeaderType>
void VtkReader<Grid,Creator,Field>::readPointsAppended (MetaType<FloatType>, MetaType<HeaderType>, std::ifstream& input)
{
  VTK_ASSERT(numberOfPoints_ > 0);
  VTK_ASSERT(dataArray_["Points"].components == 3u);
  std::vector<FloatType> point_values;
  readAppended(input, point_values, HeaderType(dataArray_["Points"].offset));
  VTK_ASSERT(point_values.size() == 3*numberOfPoints_);

  // extract points from continuous values
  GlobalCoordinate p;
  vec_points.reserve(numberOfPoints_);
  std::size_t idx = 0;
  for (std::size_t i = 0; i < numberOfPoints_; ++i) {
    for (std::size_t j = 0; j < p.size(); ++j)
      p[j] = FloatType(point_values[idx++]);
    idx += (3u - p.size());
    vec_points.push_back(p);
  }
}


template <class Grid, class Creator, class Field>
typename VtkReader<Grid,Creator,Field>::Sections
VtkReader<Grid,Creator,Field>::readCells (std::ifstream& input, std::string id)
{
  Sections sec = CELLS_DATA_ARRAY;

  VTK_ASSERT(numberOfCells_ > 0);
  if (id == "Cells.types") {
    sec = readDataArray(input, vec_types, numberOfCells_, CELLS_DATA_ARRAY, CELLS);
    VTK_ASSERT(vec_types.size() == numberOfCells_);
  } else if (id == "Cells.offsets") {
    sec = readDataArray(input, vec_offsets, numberOfCells_, CELLS_DATA_ARRAY, CELLS);
    VTK_ASSERT(vec_offsets.size() == numberOfCells_);
  } else if (id == "Cells.connectivity") {
    sec = readDataArray(input, vec_connectivity, std::size_t(-1), CELLS_DATA_ARRAY, CELLS);
  } else if (id == "Cells.global_point_ids") {
    sec = readDataArray(input, vec_point_ids, numberOfPoints_, CELLS_DATA_ARRAY, CELLS);
    VTK_ASSERT(vec_point_ids.size() == numberOfPoints_);
  }

  return sec;
}


template <class Grid, class Creator, class Field>
  template <class HeaderType>
void VtkReader<Grid,Creator,Field>::readCellsAppended (MetaType<HeaderType>, std::ifstream& input)
{
  VTK_ASSERT(numberOfCells_ > 0);
  auto types_data = dataArray_["Cells.types"];
  auto offsets_data = dataArray_["Cells.offsets"];
  auto connectivity_data = dataArray_["Cells.connectivity"];

  VTK_ASSERT(types_data.type == Vtk::DataTypes::UINT8);
  readAppended(input, vec_types, HeaderType(types_data.offset));
  VTK_ASSERT(vec_types.size() == numberOfCells_);

  if (offsets_data.type == Vtk::DataTypes::INT64)
    readAppended(input, vec_offsets, HeaderType(offsets_data.offset));
  else if (offsets_data.type == Vtk::DataTypes::INT32) {
    std::vector<std::int32_t> offsets;
    readAppended(input, offsets, HeaderType(offsets_data.offset));
    vec_offsets.resize(offsets.size());
    std::copy(offsets.begin(), offsets.end(), vec_offsets.begin());
  }
  else { DUNE_THROW(Dune::NotImplemented, "Unsupported DataType in Cell offsets."); }
  VTK_ASSERT(vec_offsets.size() == numberOfCells_);

  if (connectivity_data.type == Vtk::DataTypes::INT64)
    readAppended(input, vec_connectivity, HeaderType(connectivity_data.offset));
  else if (connectivity_data.type == Vtk::DataTypes::INT32) {
    std::vector<std::int32_t> connectivity;
    readAppended(input, connectivity, HeaderType(connectivity_data.offset));
    vec_connectivity.resize(connectivity.size());
    std::copy(connectivity.begin(), connectivity.end(), vec_connectivity.begin());
  }
  else { DUNE_THROW(Dune::NotImplemented, "Unsupported DataType in Cell connectivity."); }
  VTK_ASSERT(vec_connectivity.size() == std::size_t(vec_offsets.back()));

  if (dataArray_.count("Cells.global_point_ids") > 0) {
    auto point_id_data = dataArray_["Cells.global_point_ids"];
    VTK_ASSERT(point_id_data.type == Vtk::DataTypes::UINT64);
    readAppended(input, vec_point_ids, HeaderType(point_id_data.offset));
    VTK_ASSERT(vec_point_ids.size() == numberOfPoints_);
  }
}


// @{ implementation detail
namespace {

/**
 * Read compressed data into `buffer_in`, uncompress it and store the result in
 * the concrete-data-type `buffer`
 * \param bs     Size of the uncompressed data
 * \param cbs    Size of the compressed data
 * \param input  Stream to read from.
 **/
template <class T, class IStream>
void read_compressed_zlib (T* buffer, unsigned char* buffer_in,
                           std::uint64_t bs, std::uint64_t cbs, IStream& input)
{
#if HAVE_VTK_ZLIB
  uLongf uncompressed_space = uLongf(bs);
  uLongf compressed_space = uLongf(cbs);

  Bytef* compressed_buffer = reinterpret_cast<Bytef*>(buffer_in);
  Bytef* uncompressed_buffer = reinterpret_cast<Bytef*>(buffer);

  input.read((char*)(compressed_buffer), compressed_space);
  VTK_ASSERT(uLongf(input.gcount()) == compressed_space);

  if (uncompress(uncompressed_buffer, &uncompressed_space, compressed_buffer, compressed_space) != Z_OK) {
    std::cerr << "Zlib error while uncompressing data.\n";
    std::abort();
  }
  VTK_ASSERT(uLongf(bs) == uncompressed_space);
#else
  std::cerr << "ZLib Compression not supported. Provide the ZLIB package to CMake." << std::endl;
  std::abort();
#endif
}

template <class T, class IStream>
void read_compressed_lz4 (T* /* buffer */, unsigned char* /* buffer_in */,
                          std::uint64_t /* bs */, std::uint64_t /* cbs */, IStream& /* input */)
{
#if HAVE_VTK_LZ4
  std::cerr << "LZ4 Compression not yet implemented" << std::endl;
  std::abort();
#else
  std::cerr << "LZ4 Compression not supported. Provide the LZ4 package to CMake." << std::endl;
  std::abort();
#endif
}

template <class T, class IStream>
void read_compressed_lzma (T* /* buffer */, unsigned char* /* buffer_in */,
                           std::uint64_t /* bs */, std::uint64_t /* cbs */, IStream& /* input */)
{
#if HAVE_VTK_LZMA
  std::cerr << "LZMA Compression not yet implemented" << std::endl;
  std::abort();
#else
  std::cerr << "LZMA Compression not supported. Provide the LZMA package to CMake." << std::endl;
  std::abort();
#endif
}

}
// @}


template <class Grid, class Creator, class Field>
  template <class FloatType, class HeaderType>
void VtkReader<Grid,Creator,Field>::readAppended (std::ifstream& input, std::vector<FloatType>& values, HeaderType offset)
{
  input.seekg(offset0_ + offset);

  HeaderType size = 0;

  HeaderType num_blocks = 0;
  HeaderType block_size = 0;
  HeaderType last_block_size = 0;
  std::vector<HeaderType> cbs; // compressed block sizes

  // read total size / block-size(s)
  if (compressor_ != Vtk::CompressorTypes::NONE) {
    input.read((char*)&num_blocks, sizeof(HeaderType));
    input.read((char*)&block_size, sizeof(HeaderType));
    input.read((char*)&last_block_size, sizeof(HeaderType));

    VTK_ASSERT(block_size % sizeof(FloatType) == 0);

    // total size of the uncompressed data
    size = block_size * (num_blocks-1) + last_block_size;

    // size of the compressed blocks
    cbs.resize(num_blocks);
    input.read((char*)cbs.data(), num_blocks*sizeof(HeaderType));
  } else {
    input.read((char*)&size, sizeof(HeaderType));
  }
  VTK_ASSERT(size > 0 && (size % sizeof(FloatType)) == 0);
  values.resize(size / sizeof(FloatType));

  if (compressor_ != Vtk::CompressorTypes::NONE) {
    // upper bound for compressed block-size
    HeaderType compressed_block_size = block_size + (block_size + 999)/1000 + 12;
    // number of values in the full blocks
    std::size_t num_values = block_size / sizeof(FloatType);

    std::vector<unsigned char> buffer_in(compressed_block_size);
    for (std::size_t i = 0; i < std::size_t(num_blocks); ++i) {
      HeaderType bs = i < std::size_t(num_blocks-1) ? block_size : last_block_size;

      switch (compressor_) {
        case Vtk::CompressorTypes::ZLIB:
          read_compressed_zlib(values.data() + i*num_values, buffer_in.data(), bs, cbs[i], input);
          break;
        case Vtk::CompressorTypes::LZ4:
          read_compressed_lz4(values.data() + i*num_values, buffer_in.data(), bs, cbs[i], input);
          break;
        case Vtk::CompressorTypes::LZMA:
          read_compressed_lzma(values.data() + i*num_values, buffer_in.data(), bs, cbs[i], input);
          break;
        default:
          VTK_ASSERT_MSG(false, "Unsupported Compressor type.");
          break;
      }
    }
  } else {
    input.read((char*)(values.data()), size);
    VTK_ASSERT(input.gcount() == std::streamsize(size));
  }
}


template <class Grid, class Creator, class Field>
void VtkReader<Grid,Creator,Field>::fillGridCreator (bool insertPieces)
{
  VTK_ASSERT(vec_points.size() == numberOfPoints_);
  VTK_ASSERT(vec_types.size() == numberOfCells_);
  VTK_ASSERT(vec_offsets.size() == numberOfCells_);

  if (!vec_points.empty())
    creator_->insertVertices(vec_points, vec_point_ids);
  if (!vec_types.empty())
    creator_->insertElements(vec_types, vec_offsets, vec_connectivity);
  if (insertPieces)
    creator_->insertPieces(pieces_);
}


// Convert section into string
template <class Grid, class Creator, class Field>
std::string VtkReader<Grid,Creator,Field>::toString(Sections s) const
{
  switch (s) {
    case VTK_FILE:
      return "VTKFile";
    case UNSTRUCTURED_GRID:
      return "UnstructuredGrid";
    case PIECE:
      return "Piece";
    case POINT_DATA:
      return "PointData";
    case CELL_DATA:
      return "CellData";
    case POINTS:
      return "Points";
    case CELLS:
      return "Cells";
    case APPENDED_DATA:
      return "AppendedData";
    case PD_DATA_ARRAY:
    case CD_DATA_ARRAY:
    case POINTS_DATA_ARRAY:
    case CELLS_DATA_ARRAY:
      return "DataArray";
    default:
      return "Unknown";
  }
}


// Assume input already read the line <AppendedData ...>
template <class Grid, class Creator, class Field>
std::uint64_t VtkReader<Grid,Creator,Field>::findAppendedDataPosition (std::ifstream& input) const
{
  char c;
  while (input.get(c) && std::isblank(c)) { /*do nothing*/ }

  std::uint64_t offset = input.tellg();
  if (c != '_')
    --offset; // if char is not '_', assume it is part of the data.

  return offset;
}


template <class Grid, class Creator, class Field>
std::map<std::string, std::string> VtkReader<Grid,Creator,Field>::parseXml (std::string const& line, bool& closed)
{
  closed = false;
  std::map<std::string, std::string> attr;

  Sections sec = NO_SECTION;
  bool escape = false;

  std::string name = "";
  std::string value = "";
  for (auto c : line) {
    switch (sec) {
    case NO_SECTION:
      if (std::isalpha(c) || c == '_') {
        name.clear();
        sec = XML_NAME;
        name.push_back(c);
      } else if (c == '/') {
        closed = true;
      }
      break;
    case XML_NAME:
      if (std::isalpha(c) || c == '_')
        name.push_back(c);
      else
        sec = (c == '=' ? XML_NAME_ASSIGN : NO_SECTION);
      break;
    case XML_NAME_ASSIGN:
      value.clear();
      escape = false;
      VTK_ASSERT_MSG( c == '"', "Format error!" );
      sec = XML_VALUE;
      break;
    case XML_VALUE:
      if (c == '"' && !escape) {
        attr[name] = value;
        sec = NO_SECTION;
      } else if (c == '\\' && !escape) {
        escape = true;
      }  else {
        value.push_back(c);
        escape = false;
      }
      break;
    default:
      VTK_ASSERT_MSG(false, "Format error!");
    }
  }

  return attr;
}


template <class Grid, class Creator, class Field>
void VtkReader<Grid,Creator,Field>::clear ()
{
  vec_points.clear();
  vec_point_ids.clear();
  vec_types.clear();
  vec_offsets.clear();
  vec_connectivity.clear();
  dataArray_.clear();
  pieces_.clear();

  numberOfCells_ = 0;
  numberOfPoints_ = 0;
  offset0_ = 0;
  read_ = false;
}

} // end namespace Dune
