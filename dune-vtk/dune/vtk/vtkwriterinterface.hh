#pragma once

#include <iosfwd>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/vtk/filewriter.hh>
#include <dune/vtk/function.hh>
#include <dune/vtk/types.hh>

namespace Dune
{
  /// Interface for file writers for the Vtk XML file formats
  /**
   * \tparam GV  Model of Dune::GridView
   * \tparam DC  Model of \ref DataCollectorInterface
   **/
  template <class GV, class DC>
  class VtkWriterInterface
      : public Vtk::FileWriter
  {
    template <class> friend class TimeseriesWriter;
    template <class> friend class PvdWriter;

  public:
    using GridView = GV;
    using DataCollector = DC;

  protected:
    using VtkFunction = Dune::Vtk::Function<GridView>;
    using pos_type = typename std::ostream::pos_type;

    enum PositionTypes {
      POINT_DATA,
      CELL_DATA
    };

  public:
    /// \brief Constructor, passes the gridView to the DataCollector
    /**
     * Creates a new VtkWriterInterface for the provided GridView. Initializes a
     * DataCollector that is used to collect point coordinates, cell connectivity and
     * data values.
     *
     * This constructor assumes, that the DataCollector can be constructed from a single argument,
     * the passed gridView.
     *
     * \param gridView  Implementation of Dune::GridView
     * \param format    Format of the VTK file, either Vtk::FormatTypes::BINARY, Vtk::FormatTypes::ASCII, or Vtk::COMPRESSED
     * \param datatype  Data type of a single component of the point coordinates [Vtk::DataTypes::FLOAT32]
     * \param headertype  Integer type used in binary data headers [Vtk::DataTypes::UINT32]
     **/
    VtkWriterInterface (GridView const& gridView,
                        Vtk::FormatTypes format = Vtk::FormatTypes::BINARY,
                        Vtk::DataTypes datatype = Vtk::DataTypes::FLOAT32,
                        Vtk::DataTypes headertype = Vtk::DataTypes::UINT32)
      : VtkWriterInterface(std::make_shared<DataCollector>(gridView), format, datatype, headertype)
    {}

    /// \brief Constructor, wraps the passed DataCollector in a non-destroying shared_ptr
    VtkWriterInterface (DataCollector& dataCollector,
                        Vtk::FormatTypes format = Vtk::FormatTypes::BINARY,
                        Vtk::DataTypes datatype = Vtk::DataTypes::FLOAT32,
                        Vtk::DataTypes headertype = Vtk::DataTypes::UINT32)
      : VtkWriterInterface(stackobject_to_shared_ptr(dataCollector), format, datatype, headertype)
    {}

    /// \brief Constructor, stores the passed DataCollector
    VtkWriterInterface (std::shared_ptr<DataCollector> dataCollector,
                        Vtk::FormatTypes format = Vtk::FormatTypes::BINARY,
                        Vtk::DataTypes datatype = Vtk::DataTypes::FLOAT32,
                        Vtk::DataTypes headertype = Vtk::DataTypes::UINT32)
      : dataCollector_(std::move(dataCollector))
    {
      setFormat(format);
      setDatatype(datatype);
      setHeadertype(headertype);
    }


    /// \brief Write the attached data to the file
    /**
     * \param fn   Filename of the VTK file. May contain a directory and any file extension.
     * \param dir  The optional parameter specifies the directory of the partition files for parallel writes.
     *
     * \returns File name that is actually written.
     **/
    virtual std::string write (std::string const& fn, std::optional<std::string> dir = {}) const override;

    /// \brief Attach point data to the writer
    /**
     * Attach a global function to the writer that will be evaluated at grid points
     * (vertices and higher order points). The global function must be
     * assignable to the function wrapper \ref Vtk::Function. Additional argument
     * for output datatype and number of components can be passed. See \ref Vtk::Function
     * Constructor for possible arguments.
     *
     * \param fct     A GridFunction, LocalFunction, or Dune::VTKFunction
     * \param args... Additional arguments, like `name`, `numComponents`, `dataType` or `Vtk::FieldInfo`
     **/
    template <class Function, class... Args>
    VtkWriterInterface& addPointData (Function&& fct, Args&&... args)
    {
      pointData_.emplace_back(std::forward<Function>(fct), std::forward<Args>(args)...,
                              datatype_, Vtk::RangeTypes::AUTO);
      return *this;
    }

    /// \brief Attach cell data to the writer
    /**
     * Attach a global function to the writer that will be evaluated at cell centers.
     * The global function must be assignable to the function wrapper \ref Vtk::Function.
     * Additional argument for output datatype and number of components can be passed.
     * See \ref Vtk::Function Constructor for possible arguments.
     *
     * \param fct     A GridFunction, LocalFunction, or Dune::VTKFunction
     * \param args... Additional arguments, like `name`, `numComponents`, `dataType` or `Vtk::FieldInfo`
     **/
    template <class Function, class... Args>
    VtkWriterInterface& addCellData (Function&& fct, Args&&... args)
    {
      cellData_.emplace_back(std::forward<Function>(fct), std::forward<Args>(args)...,
                             datatype_, Vtk::RangeTypes::AUTO);
      return *this;
    }


    // Sets the VTK file format
    void setFormat (Vtk::FormatTypes format)
    {
      format_ = format;

      if (format_ == Vtk::FormatTypes::COMPRESSED) {
#if HAVE_VTK_ZLIB
        compressor_ = Vtk::CompressorTypes::ZLIB;
#else
        std::cout << "Dune is compiled without compression. Falling back to BINARY VTK output!\n";
        format_ = Vtk::FormatTypes::BINARY;
#endif
      } else {
        compressor_ = Vtk::CompressorTypes::NONE;
      }

    }

    /// Sets the global datatype used for coordinates and other global float values
    void setDatatype (Vtk::DataTypes datatype)
    {
      datatype_ = datatype;
    }

    /// Sets the integer type used in binary data headers
    void setHeadertype (Vtk::DataTypes datatype)
    {
      headertype_ = datatype;
    }

    /// Sets the compressor type used in binary data headers, Additionally a compression
    /// level can be passed with level = -1 means: default compression level. Level must be in [0-9]
    void setCompressor (Vtk::CompressorTypes compressor, int level = -1)
    {
      compressor_ = compressor;
      compression_level = level;
      VTK_ASSERT(level >= -1 && level <= 9);

      if (compressor_ != Vtk::CompressorTypes::NONE)
        format_ = Vtk::FormatTypes::COMPRESSED;
    }

  private:
    /// Write a serial VTK file in Unstructured format
    virtual void writeSerialFile (std::ofstream& out) const = 0;

    /// Write a parallel VTK file `pfilename.pvtx` in XML format,
    /// with `size` the number of pieces and serial files given by `pfilename_p[i].vtu`
    /// for [i] in [0,...,size).
    virtual void writeParallelFile (std::ofstream& out, std::string const& pfilename, int size) const = 0;

    /// Return the file extension of the serial file (not including the dot)
    virtual std::string fileExtension () const = 0;

    /// Write points and cells in raw/compressed format to output stream
    virtual void writeGridAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const = 0;

  protected:
    // Write the point or cell values given by the grid function `fct` to the
    // output stream `out`. In case of binary format, append the streampos of XML
    // attributes "offset" to the vector `offsets`.
    void writeData (std::ofstream& out,
                    std::vector<pos_type>& offsets,
                    VtkFunction const& fct,
                    PositionTypes type,
                    std::optional<std::size_t> timestep = {}) const;

    // Write point-data and cell-data in raw/compressed format to output stream
    void writeDataAppended (std::ofstream& out, std::vector<std::uint64_t>& blocks) const;

    // Write the coordinates of the vertices to the output stream `out`. In case
    // of binary format, appends the streampos of XML attributes "offset" to the
    // vector `offsets`.
    void writePoints (std::ofstream& out,
                      std::vector<pos_type>& offsets,
                      std::optional<std::size_t> timestep = {}) const;

    // Write Appended section and fillin offset values to XML attributes
    void writeAppended (std::ofstream& out, std::vector<pos_type> const& offsets) const;

    // Write the `values` in blocks (possibly compressed) to the output
    // stream `out`. Return the written block size.
    template <class HeaderType, class FloatType>
    std::uint64_t writeValuesAppended (std::ofstream& out, std::vector<FloatType> const& values) const;

    // Write the `values` in a space and newline separated list of ascii representations.
    // The precision is controlled by the datatype and numerical_limits::digits10.
    template <class T>
    void writeValuesAscii (std::ofstream& out, std::vector<T> const& values) const;

    // Write the XML file header of a VTK file `<VTKFile ...>`
    void writeHeader (std::ofstream& out, std::string const& type) const;

    /// Return PointData/CellData attributes for the name of the first scalar/vector/tensor DataArray
    std::string getNames (std::vector<VtkFunction> const& data) const;

    // Returns endianness
    std::string getEndian () const
    {
      short i = 1;
      return (reinterpret_cast<char*>(&i)[1] == 1 ? "BigEndian" : "LittleEndian");
    }

    // provide accessor to \ref fileExtension virtual method
    std::string getFileExtension () const
    {
      return fileExtension();
    }

    // Returns the VTK file format initialized in the constructor
    Vtk::FormatTypes getFormat () const
    {
      return format_;
    }

    // Returns the global datatype used for coordinates and other global float values
    Vtk::DataTypes getDatatype () const
    {
      return datatype_;
    }

    // Return the global MPI communicator.
    auto comm () const
    {
      return MPIHelper::getCollectiveCommunication();
    }

  protected:
    std::shared_ptr<DataCollector> dataCollector_;

    Vtk::FormatTypes format_;
    Vtk::DataTypes datatype_;
    Vtk::DataTypes headertype_;
    Vtk::CompressorTypes compressor_ = Vtk::CompressorTypes::NONE;

    // attached data
    std::vector<VtkFunction> pointData_;
    std::vector<VtkFunction> cellData_;

    std::size_t const block_size = 1024*32;
    int compression_level = -1; // in [0,9], -1 ... use default value
  };


  template <class Writer>
  struct IsVtkWriter
  {
    template <class GV, class DC>
    static std::uint16_t test(VtkWriterInterface<GV,DC> const&);
    static std::uint8_t  test(...); // fall-back overload

    static constexpr bool value = sizeof(test(std::declval<Writer>())) > sizeof(std::uint8_t);
  };

} // end namespace Dune

#include "vtkwriterinterface.impl.hh"
