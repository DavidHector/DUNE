#pragma once

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#if HAVE_VTK_ZLIB
#include <zlib.h>
#endif

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/vtk/utility/enum.hh>
#include <dune/vtk/utility/errors.hh>
#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/string.hh>

namespace Dune {

template <class GV, class DC>
std::string VtkWriterInterface<GV,DC>
  ::write (std::string const& fn, std::optional<std::string> dir) const
{
  dataCollector_->update();

  auto p = Vtk::Path(fn);
  auto name = p.stem();
  p.removeFilename();

  Vtk::Path fn_dir = p;
  Vtk::Path data_dir = dir ? Vtk::Path(*dir) : fn_dir;
  Vtk::Path rel_dir = Vtk::relative(data_dir, fn_dir);

  std::string serial_fn = data_dir.string() + '/' + name.string();
  std::string parallel_fn = fn_dir.string() + '/' + name.string();
  std::string rel_fn = rel_dir.string() + '/' + name.string();

  if (comm().size() > 1)
    serial_fn += "_p" + std::to_string(comm().rank());

  std::string outputFilename;

  { // write serial file
    outputFilename = serial_fn + "." + fileExtension();
    std::ofstream serial_out(outputFilename, std::ios_base::ate | std::ios::binary);
    assert(serial_out.is_open());

    serial_out.imbue(std::locale::classic());
    serial_out << std::setprecision(datatype_ == Vtk::DataTypes::FLOAT32
      ? std::numeric_limits<float>::max_digits10
      : std::numeric_limits<double>::max_digits10);

    writeSerialFile(serial_out);
  }

  if (comm().size() > 1 && comm().rank() == 0) {
    // write parallel filee
    outputFilename = parallel_fn + ".p" + fileExtension();
    std::ofstream parallel_out(outputFilename, std::ios_base::ate | std::ios::binary);
    assert(parallel_out.is_open());

    parallel_out.imbue(std::locale::classic());
    parallel_out << std::setprecision(datatype_ == Vtk::DataTypes::FLOAT32
      ? std::numeric_limits<float>::max_digits10
      : std::numeric_limits<double>::max_digits10);

    writeParallelFile(parallel_out, rel_fn, comm().size());
  }

  return outputFilename;
}


template <class GV, class DC>
void VtkWriterInterface<GV,DC>
  ::writeData (std::ofstream& out, std::vector<pos_type>& offsets,
               VtkFunction const& fct, PositionTypes type,
               std::optional<std::size_t> timestep) const
{
  out << "<DataArray Name=\"" << fct.name() << "\" type=\"" << to_string(fct.dataType()) << "\""
      << " NumberOfComponents=\"" << fct.numComponents() << "\" format=\"" << (format_ == Vtk::FormatTypes::ASCII ? "ascii\"" : "appended\"");
  if (timestep)
    out << " TimeStep=\"" << *timestep << "\"";

  if (format_ == Vtk::FormatTypes::ASCII) {
    out << ">\n";
    if (type == POINT_DATA)
      writeValuesAscii(out, dataCollector_->template pointData<double>(fct));
    else
      writeValuesAscii(out, dataCollector_->template cellData<double>(fct));
    out << "</DataArray>\n";
  } else {
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";
  }
}


template <class GV, class DC>
void VtkWriterInterface<GV,DC>
  ::writePoints (std::ofstream& out, std::vector<pos_type>& offsets,
                std::optional<std::size_t> timestep) const
{
  out << "<DataArray type=\"" << to_string(datatype_) << "\""
      << " NumberOfComponents=\"3\" format=\"" << (format_ == Vtk::FormatTypes::ASCII ? "ascii\"" : "appended\"");
  if (timestep)
    out << " TimeStep=\"" << *timestep << "\"";

  if (format_ == Vtk::FormatTypes::ASCII) {
    out << ">\n";
    writeValuesAscii(out, dataCollector_->template points<double>());
    out << "</DataArray>\n";
  } else {
    out << " offset=";
    offsets.push_back(out.tellp());
    out << std::string(std::numeric_limits<std::uint64_t>::digits10 + 2, ' ');
    out << "/>\n";
  }
}


template <class GV, class DC>
void VtkWriterInterface<GV,DC>
  ::writeDataAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const
{
  for (auto const& v : pointData_) {
    Vtk::mapDataTypes<std::is_floating_point, std::is_integral>(v.dataType(), headertype_,
    [&](auto f, auto h) {
      using F = typename decltype(f)::type;
      using H = typename decltype(h)::type;
      blocks.push_back(this->template writeValuesAppended<H>(out, dataCollector_->template pointData<F>(v)));
    });
  }

  for (auto const& v : cellData_) {
    Vtk::mapDataTypes<std::is_floating_point, std::is_integral>(v.dataType(), headertype_,
    [&](auto f, auto h) {
      using F = typename decltype(f)::type;
      using H = typename decltype(h)::type;
      blocks.push_back(this->template writeValuesAppended<H>(out, dataCollector_->template cellData<F>(v)));
    });
  }
}


template <class GV, class DC>
void VtkWriterInterface<GV,DC>
  ::writeAppended (std::ofstream& out, std::vector<pos_type> const& offsets) const
{
  if (is_a(format_, Vtk::FormatTypes::APPENDED)) {
    out << "<AppendedData encoding=\"raw\">\n_";
    std::vector<std::uint64_t> blocks;
    writeGridAppended(out, blocks);
    writeDataAppended(out, blocks);
    out << "</AppendedData>\n";
    pos_type appended_pos = out.tellp();

    pos_type offset = 0;
    for (std::size_t i = 0; i < offsets.size(); ++i) {
      out.seekp(offsets[i]);
      out << '"' << offset << '"';
      offset += pos_type(blocks[i]);
    }

    out.seekp(appended_pos);
  }
}


namespace Impl {

  template <class T, std::enable_if_t<(sizeof(T)>1), int> = 0>
  inline T const& printable (T const& t) { return t; }

  inline std::int16_t printable (std::int8_t c) { return std::int16_t(c); }
  inline std::uint16_t printable (std::uint8_t c) { return std::uint16_t(c); }

} // end namespace Impl


template <class GV, class DC>
  template <class T>
void VtkWriterInterface<GV,DC>
  ::writeValuesAscii (std::ofstream& out, std::vector<T> const& values) const
{
  assert(is_a(format_, Vtk::FormatTypes::ASCII) && "Function should by called only in ascii mode!\n");
  std::size_t i = 0;
  for (auto const& v : values)
    out << Impl::printable(v) << (++i % 6 != 0 ? ' ' : '\n');
  if (i % 6 != 0)
    out << '\n';
}

template <class GV, class DC>
void VtkWriterInterface<GV,DC>
  ::writeHeader (std::ofstream& out, std::string const& type) const
{
  out << "<VTKFile"
      << " type=\"" << type << "\""
      << " version=\"1.0\""
      << " header_type=\"" << to_string(headertype_) << "\"";
  if (format_ != Vtk::FormatTypes::ASCII)
    out << " byte_order=\"" << getEndian() << "\"";
  if (compressor_ != Vtk::CompressorTypes::NONE)
    out << " compressor=\"" << to_string(compressor_) << "\"";
  out << ">\n";
}


namespace Impl {

template <class T>
std::uint64_t writeValuesToBuffer (std::size_t max_num_values, unsigned char* buffer,
                                   std::vector<T> const& vec, std::size_t shift)
{
  std::size_t num_values = std::min(max_num_values, vec.size()-shift);
  std::uint64_t bs = num_values*sizeof(T);
  std::memcpy(buffer, (unsigned char*)(vec.data()+shift), std::size_t(bs));
  return bs;
}

inline std::uint64_t compressBuffer_zlib (unsigned char const* buffer, unsigned char* buffer_out,
                                          std::uint64_t bs, std::uint64_t cbs, int level)
{
#if HAVE_VTK_ZLIB
  uLongf uncompressed_space = uLongf(bs);
  uLongf compressed_space = uLongf(cbs);

  Bytef* out = reinterpret_cast<Bytef*>(buffer_out);
  Bytef const* in = reinterpret_cast<Bytef const*>(buffer);

  if (compress2(out, &compressed_space, in, uncompressed_space, level) != Z_OK) {
    std::cerr << "Zlib error while compressing data.\n";
    std::abort();
  }

  return compressed_space;
#else
  std::cerr << "Can not call writeCompressed without compression enabled!\n";
  std::abort();
  return 0;
#endif
}

inline std::uint64_t compressBuffer_lz4 (unsigned char const* /* buffer */, unsigned char* /* buffer_out */,
                                         std::uint64_t /* bs */, std::uint64_t /* cbs */, int /* level */)
{
#if HAVE_VTK_LZ4
  std::cerr << "LZ4 Compression not yet implemented" << std::endl;
  std::abort();
  return 0;
#else
  std::cerr << "LZ4 Compression not supported. Provide the LZ4 package to CMake." << std::endl;
  std::abort();
  return 0;
#endif
}

inline std::uint64_t compressBuffer_lzma (unsigned char const* /* buffer */, unsigned char* /* buffer_out */,
                                          std::uint64_t /* bs */, std::uint64_t /* cbs */, int /* level */)
{
#if HAVE_VTK_LZMA
  std::cerr << "LZMA Compression not yet implemented" << std::endl;
  std::abort();
  return 0;
#else
  std::cerr << "LZMA Compression not supported. Provide the LZMA package to CMake." << std::endl;
  std::abort();
  return 0;
#endif
}

} // end namespace Impl

template <class GV, class DC>
  template <class HeaderType, class FloatType>
std::uint64_t VtkWriterInterface<GV,DC>
  ::writeValuesAppended (std::ofstream& out, std::vector<FloatType> const& values) const
{
  assert(is_a(format_, Vtk::FormatTypes::APPENDED) && "Function should by called only in appended mode!\n");
  pos_type begin_pos = out.tellp();

  HeaderType size = values.size() * sizeof(FloatType);

  HeaderType num_full_blocks = size / block_size;
  HeaderType last_block_size = size % block_size;
  HeaderType num_blocks = num_full_blocks + (last_block_size > 0 ? 1 : 0);

  // write block-size(s)
  HeaderType zero = 0;
  if (compressor_ != Vtk::CompressorTypes::NONE) {
    out.write((char*)&num_blocks, sizeof(HeaderType));
    out.write((char*)&block_size, sizeof(HeaderType));
    out.write((char*)&last_block_size, sizeof(HeaderType));
    for (HeaderType i = 0; i < num_blocks; ++i)
      out.write((char*)&zero, sizeof(HeaderType));
  } else {
    out.write((char*)&size, sizeof(HeaderType));
  }

  HeaderType compressed_block_size = block_size + (block_size + 999)/1000 + 12;
  std::vector<unsigned char> buffer(block_size);
  std::vector<unsigned char> buffer_out;
  if (compressor_ != Vtk::CompressorTypes::NONE)
    buffer_out.resize(std::size_t(compressed_block_size));

  std::size_t num_values = block_size / sizeof(FloatType);

  std::vector<HeaderType> cbs(std::size_t(num_blocks), 0); // compressed block sizes
  for (std::size_t i = 0; i < std::size_t(num_blocks); ++i) {
    HeaderType bs = Impl::writeValuesToBuffer<FloatType>(num_values, buffer.data(), values, i*num_values);

    switch (compressor_) {
      case Vtk::CompressorTypes::NONE:
        out.write((char*)buffer.data(), bs);
        break;
      case Vtk::CompressorTypes::ZLIB:
        cbs[i] = Impl::compressBuffer_zlib(buffer.data(), buffer_out.data(), bs, compressed_block_size, compression_level);
        out.write((char*)buffer_out.data(), cbs[i]);
        break;
      case Vtk::CompressorTypes::LZ4:
        cbs[i] = Impl::compressBuffer_zlib(buffer.data(), buffer_out.data(), bs, compressed_block_size, compression_level);
        out.write((char*)buffer_out.data(), cbs[i]);
        break;
      case Vtk::CompressorTypes::LZMA:
        cbs[i] = Impl::compressBuffer_zlib(buffer.data(), buffer_out.data(), bs, compressed_block_size, compression_level);
        out.write((char*)buffer_out.data(), cbs[i]);
        break;
    }
  }

  pos_type end_pos = out.tellp();
  if (compressor_ != Vtk::CompressorTypes::NONE) {
    out.seekp(begin_pos + std::streamoff(3*sizeof(HeaderType)));
    out.write((char*)cbs.data(), std::streamsize(num_blocks*sizeof(HeaderType)));
    out.seekp(end_pos);
  }

  return std::uint64_t(end_pos - begin_pos);
}


template <class GV, class DC>
std::string VtkWriterInterface<GV,DC>
  ::getNames (std::vector<VtkFunction> const& data) const
{
  auto scalar = std::find_if(data.begin(), data.end(), [](auto const& v) { return v.rangeType() == Vtk::RangeTypes::SCALAR; });
  auto vector = std::find_if(data.begin(), data.end(), [](auto const& v) { return v.rangeType() == Vtk::RangeTypes::VECTOR; });
  auto tensor = std::find_if(data.begin(), data.end(), [](auto const& v) { return v.rangeType() == Vtk::RangeTypes::TENSOR; });
  return (scalar != data.end() ? " Scalars=\"" + scalar->name() + "\"" : "")
       + (vector != data.end() ? " Vectors=\"" + vector->name() + "\"" : "")
       + (tensor != data.end() ? " Tensors=\"" + tensor->name() + "\"" : "");
}

} // end namespace Dune
