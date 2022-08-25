#pragma once

#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>

#if HAVE_VTK_ZLIB
#include <zlib.h>
#endif

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/vtk/utility/enum.hh>
#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/string.hh>

namespace Dune {

template <class W>
VtkTimeseriesWriter<W>::~VtkTimeseriesWriter ()
{
  if (initialized_) {
    [[maybe_unused]] int ec = std::remove(filenameMesh_.c_str());
    assert(ec == 0);
    for (auto const& timestep : timesteps_) {
      ec = std::remove(timestep.second.c_str());
      assert(ec == 0);
    }
  }
  std::remove(tmpDir_.string().c_str());
}


template <class W>
void VtkTimeseriesWriter<W>
  ::writeTimestep (double time, std::string const& fn, std::optional<std::string> tmpDir, bool writeCollection) const
{
  auto name = Vtk::Path(fn).stem();
  auto tmpBase = tmpDir ? Vtk::Path(*tmpDir) : tmpDir_;
  auto tmp = tmpBase;
  tmp /= name.string();

  vtkWriter_.dataCollector_->update();

  std::string filenameBase = tmp.string();

  if (vtkWriter_.comm().size() > 1)
    filenameBase = tmp.string() + "_p" + std::to_string(vtkWriter_.comm().rank());

  if (!initialized_) {
    Vtk::createDirectories(tmpBase);

    // write points and cells only once
    filenameMesh_ = filenameBase + ".mesh.vtkdata";
    std::ofstream out(filenameMesh_, std::ios_base::ate | std::ios::binary);
    vtkWriter_.writeGridAppended(out, blocks_);
    initialized_ = true;
  }

  std::string filenameData = filenameBase + "_t" + std::to_string(timesteps_.size()) + ".vtkdata";
  std::ofstream out(filenameData, std::ios_base::ate | std::ios::binary);
  vtkWriter_.writeDataAppended(out, blocks_);
  timesteps_.emplace_back(time, filenameData);

  if (writeCollection)
    write(fn);
}


template <class W>
std::string VtkTimeseriesWriter<W>
  ::write (std::string const& fn, std::optional<std::string> dir) const
{
  assert( initialized_ );

  auto p = Vtk::Path(fn);
  auto name = p.stem();
  p.removeFilename();

  Vtk::Path fn_dir = p;
  Vtk::Path data_dir = dir ? Vtk::Path(*dir) : fn_dir;
  Vtk::Path rel_dir = Vtk::relative(data_dir, fn_dir);

  std::string serial_fn = fn_dir.string() + '/' + name.string() + "_ts";
  std::string parallel_fn = data_dir.string() + '/' + name.string() + "_ts";
  std::string rel_fn = rel_dir.string() + '/' + name.string() + "_ts";

  int commRank = vtkWriter_.comm().rank();
  int commSize = vtkWriter_.comm().size();
  if (commSize > 1)
    serial_fn += "_p" + std::to_string(commRank);

  std::string outputFilename;

  { // write serial file
    outputFilename = serial_fn + "." + vtkWriter_.getFileExtension();
    std::ofstream serial_out(outputFilename, std::ios_base::ate | std::ios::binary);
    assert(serial_out.is_open());

    serial_out.imbue(std::locale::classic());
    serial_out << std::setprecision(vtkWriter_.getDatatype() == Vtk::DataTypes::FLOAT32
      ? std::numeric_limits<float>::digits10+2
      : std::numeric_limits<double>::digits10+2);

    vtkWriter_.writeTimeseriesSerialFile(serial_out, filenameMesh_, timesteps_, blocks_);
  }

  if (commSize > 1 && commRank == 0) {
    // write parallel file
    outputFilename = parallel_fn + ".p" + vtkWriter_.getFileExtension();
    std::ofstream parallel_out(outputFilename, std::ios_base::ate | std::ios::binary);
    assert(parallel_out.is_open());

    parallel_out.imbue(std::locale::classic());
    parallel_out << std::setprecision(vtkWriter_.getDatatype() == Vtk::DataTypes::FLOAT32
      ? std::numeric_limits<float>::digits10+2
      : std::numeric_limits<double>::digits10+2);

    vtkWriter_.writeTimeseriesParallelFile(parallel_out, rel_fn, commSize, timesteps_);
  }

  return outputFilename;
}

} // end namespace Dune
