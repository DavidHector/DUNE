#pragma once

#include <iomanip>
#include <limits>

#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/string.hh>

namespace Dune {

template <class W>
void PvdWriter<W>
  ::writeTimestep (double time, std::string const& fn, std::optional<std::string> dir, bool writeCollection) const
{
  auto p = Vtk::Path(fn);
  auto name = p.stem();
  p.removeFilename();

  Vtk::Path fn_dir = p;
  Vtk::Path data_dir = dir ? Vtk::Path(*dir) : fn_dir;
  Vtk::Path rel_dir = Vtk::relative(data_dir, fn_dir);

  std::string pvd_fn = fn_dir.string() + '/' + name.string();
  std::string seq_fn = data_dir.string() + '/' + name.string() + "_t" + std::to_string(timesteps_.size());
  std::string rel_fn = rel_dir.string() + '/' + name.string() + "_t" + std::to_string(timesteps_.size());

  std::string ext = "." + vtkWriter_.getFileExtension();

  int commRank = vtkWriter_.comm().rank();
  int commSize = vtkWriter_.comm().size();
  if (commSize > 1)
    ext = ".p" + vtkWriter_.getFileExtension();

  timesteps_.emplace_back(time, rel_fn + ext);
  vtkWriter_.write(seq_fn + ext);

  if (commRank == 0 && writeCollection) {
    std::ofstream out(pvd_fn + ".pvd", std::ios_base::ate | std::ios::binary);
    assert(out.is_open());

    out.imbue(std::locale::classic());
    out << std::setprecision(datatype_ == Vtk::DataTypes::FLOAT32
      ? std::numeric_limits<float>::max_digits10
      : std::numeric_limits<double>::max_digits10);

    writeFile(out);
  }
}


template <class W>
std::string PvdWriter<W>
  ::write (std::string const& fn, std::optional<std::string> /*dir*/) const
{
  auto p = Vtk::Path(fn);
  auto name = p.stem();
  p.removeFilename();
  p /= name.string();

  std::string outputFilename;

  int commRank = vtkWriter_.comm().rank();
  if (commRank == 0) {
    outputFilename = p.string() + ".pvd";
    std::ofstream out(outputFilename, std::ios_base::ate | std::ios::binary);
    assert(out.is_open());

    out.imbue(std::locale::classic());
    out << std::setprecision(datatype_ == Vtk::DataTypes::FLOAT32
      ? std::numeric_limits<float>::max_digits10
      : std::numeric_limits<double>::max_digits10);

    writeFile(out);
  }

  return outputFilename;
}


template <class W>
void PvdWriter<W>
  ::writeFile (std::ofstream& out) const
{
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile"
      << " type=\"Collection\""
      << " version=\"0.1\""
      << (format_ != Vtk::FormatTypes::ASCII ? " byte_order=\"" + vtkWriter_.getEndian() + "\"" : "")
      << (format_ == Vtk::FormatTypes::COMPRESSED ? " compressor=\"vtkZLibDataCompressor\"" : "")
      << ">\n";

  out << "<Collection>\n";

  // Write all timesteps
  for (auto const& timestep : timesteps_) {
    out << "<DataSet"
        << " timestep=\"" << timestep.first << "\""
        << " part=\"0\""
        << " file=\"" << timestep.second << "\""
        << " />\n";
  }

  out << "</Collection>\n";
  out << "</VTKFile>";
}

} // end namespace Dune
