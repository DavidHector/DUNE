#pragma once

#include <iosfwd>
#include <optional>
#include <string>
#include <vector>
#include <tuple>

#include <dune/vtk/types.hh>
#include <dune/vtk/filewriter.hh>

namespace Dune
{
  /// File-Writer for ParaView .pvd files
  template <class VtkWriter>
  class PvdWriter
      : public Vtk::FileWriter
  {
    using Self = PvdWriter;

  public:
    /// Constructor, creates a VtkWriter with constructor arguments forwarded
    template <class... Args,
      disableCopyMove<Self,Args...> = 0>
    explicit PvdWriter (Args&&... args)
      : vtkWriter_{std::forward<Args>(args)...}
    {
      format_ = vtkWriter_.getFormat();
      datatype_ = vtkWriter_.getDatatype();
    }

    /// \brief Write the attached data to the file
    /**
     * Create timestep files for the data associated to the current timestep `time`.
     *
     * \param time  The time value of the written data
     * \param fn    Filename of the PVD file to write to. The base part is used to
     *              create filenames for the timestep files that are stored in \ref timesteps_.
     *              May contain directory and any filename extension.
     * \param dir   Specifies where to write the timestep files.
     * \param writeCollection  Create a collection .pvd file directly
     **/
    void writeTimestep (double time, std::string const& fn, std::optional<std::string> dir = {},
                        bool writeCollection = true) const;

    /// \brief Writes collection of timesteps to .pvd file.
    // NOTE: requires an aforegoing call to \ref writeTimestep
    /**
     * \param fn  The filename of the PVD file. May contain directory and any filename extension.
     * \param dir (Ignored) Timestep files are already written and their filenames are
     *            stored in \ref timesteps_.
     **/
    virtual std::string write (std::string const& fn, std::optional<std::string> dir = {}) const override;

    /// Attach point data to the writer, \see Vtk::Function for possible arguments
    template <class Function, class... Args>
    PvdWriter& addPointData (Function const& fct, Args&&... args)
    {
      vtkWriter_.addPointData(fct, std::forward<Args>(args)...);
      return *this;
    }

    /// Attach cell data to the writer, \see Vtk::Function for possible arguments
    template <class Function, class... Args>
    PvdWriter& addCellData (Function const& fct, Args&&... args)
    {
      vtkWriter_.addCellData(fct, std::forward<Args>(args)...);
      return *this;
    }

  protected:
    /// Write a series of vtk files in a .pvd ParaView Data file
    void writeFile (std::ofstream& out) const;

  protected:
    VtkWriter vtkWriter_;
    Vtk::FormatTypes format_;
    Vtk::DataTypes datatype_;

    mutable std::vector<std::pair<double, std::string>> timesteps_;
  };

} // end namespace Dune

#include "pvdwriter.impl.hh"
