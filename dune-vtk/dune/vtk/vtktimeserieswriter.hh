#pragma once

#include <cassert>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include <dune/vtk/filewriter.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/utility/filesystem.hh>
#include <dune/vtk/utility/uid.hh>

namespace Dune
{
  /// File-Writer for Vtk timeseries .vtu files
  /**
   * \tparam VtkWriter  Type of a FileWriter derived from \ref VtkWriterInterface that
   *                    additionally supports writeTimeseriesSerialFile() and writeTimeseriesParallelFile(),
   *                    e.g. \ref VtkUnstructuredGridWriter.
   **/
  template <class VtkWriter>
  class VtkTimeseriesWriter
      : public Vtk::FileWriter
  {
  protected:
    using Self = VtkTimeseriesWriter;
    using pos_type = typename std::ostream::pos_type;

    template <class W>
    using HasWriteTimeseriesSerialFile = decltype(&W::writeTimeseriesSerialFile);
    // static_assert(Std::is_detected<HasWriteTimeseriesSerialFile, VtkWriter>::value, "");

    template <class W>
    using HasWriteTimeseriesParallelFile = decltype(&W::writeTimeseriesParallelFile);
    // static_assert(Std::is_detected<HasWriteTimeseriesParallelFile, VtkWriter>::value, "");

  public:
    /// Constructor, stores the gridView
    template <class... Args,
      disableCopyMove<Self, Args...> = 0>
    VtkTimeseriesWriter (Args&&... args)
      : vtkWriter_{std::forward<Args>(args)...}
    {
      assert(vtkWriter_.format_ != Vtk::FormatTypes::ASCII && "Timeseries writer requires APPENDED mode");
      std::srand(std::time(nullptr));

      // put temporary file to /tmp directory
      tmpDir_ = Vtk::Path("/tmp/vtktimeserieswriter_" + Vtk::uid(10) + "/");
    }

    /// Remove all written intermediate files and remove temporary directory
    ~VtkTimeseriesWriter ();

    /// Write the attached data to the file
    /**
     * Create intermediate files for the data associated to the current timestep `time`.
     *
     * \param time  The time value of the written data
     * \param fn  Filename of the file to write to. Only the base part
     *            (without dir and extentsion) is used to write the intermediate
     *            file into a tmp directory.
     * \param tmpDir  If the directory is given, it is used as tmp dir, otherwise /tmp.
     * \param writeCollection  Create a timeseries file directly
     **/
    void writeTimestep (double time, std::string const& fn,
                        std::optional<std::string> tmpDir = {},
                        bool writeCollection = true) const;

    /// Writes all timesteps to single timeseries file.
    // NOTE: requires an aforegoing call to \ref writeTimestep
    /**
     * Create a timeseries file with all timesteps written by \ref writeTimestep.
     *
     * \param fn   Filename of the Timeseries file. May contain a directory and any file extension.
     * \param dir  The optional parameter specifies the directory of the partition files.
     *
     * \returns File name of the actual written file
     **/
    virtual std::string write (std::string const& fn, std::optional<std::string> dir = {}) const override;

    /// Attach point data to the writer, \see VtkFunction for possible arguments
    template <class Function, class... Args>
    VtkTimeseriesWriter& addPointData (Function const& fct, Args&&... args)
    {
      vtkWriter_.addPointData(fct, std::forward<Args>(args)...);
      return *this;
    }

    /// Attach cell data to the writer, \see Vtk::Function for possible arguments
    template <class Function, class... Args>
    VtkTimeseriesWriter& addCellData (Function const& fct, Args&&... args)
    {
      vtkWriter_.addCellData(fct, std::forward<Args>(args)...);
      return *this;
    }

  protected:
    VtkWriter vtkWriter_;
    Vtk::Path tmpDir_;

    mutable bool initialized_ = false;

    // block size of attached data
    mutable std::vector<std::uint64_t> blocks_;

    mutable std::string filenameMesh_;
    mutable std::vector<std::pair<double, std::string>> timesteps_;
  };

} // end namespace Dune

#include "vtktimeserieswriter.impl.hh"
