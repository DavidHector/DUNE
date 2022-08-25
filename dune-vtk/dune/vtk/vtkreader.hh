#pragma once

#include <iosfwd>
#include <map>
#include <memory>
#include <vector>

#include <dune/common/shared_ptr.hh>
#include <dune/common/typelist.hh>
#include <dune/common/typeutilities.hh>

#include <dune/vtk/filereader.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/utility/errors.hh>

// default GridCreator
#include <dune/vtk/gridcreators/continuousgridcreator.hh>
#include <dune/vtk/gridfunctions/common.hh>

namespace Dune
{
  /// File-Reader for Vtk unstructured .vtu files
  /**
   * Reads .vtu files and constructs a grid from the cells stored in the file
   * Additionally, stored data can be read.
   *
   * NOTE: Assumption on the file structure: Each XML tag must be on a separate line.
   *
   * \tparam Grid       The type of the grid to construct.
   * \tparam GC         GridCreator policy type to control what to pass to a grid factory with
   *                    data given from the file. [ContinuousGridCreator]
   * \tparam FieldType  Type of the components of the data to extract from the file [default: double]
   **/
  template <class Grid, class GC = Vtk::ContinuousGridCreator<Grid>, class FieldType = double>
  class VtkReader
      : public Vtk::FileReader<Grid, VtkReader<Grid, GC>>
  {
    // Sections visited during the xml parsing
    enum Sections {
      NO_SECTION = 0, VTK_FILE, UNSTRUCTURED_GRID, PIECE, POINT_DATA, PD_DATA_ARRAY, CELL_DATA, CD_DATA_ARRAY,
      POINTS, POINTS_DATA_ARRAY, CELLS, CELLS_DATA_ARRAY, APPENDED_DATA, XML_NAME, XML_NAME_ASSIGN, XML_VALUE
    };

    // Type storing information about read data
    struct DataArrayAttributes
    {
      std::string name;
      Vtk::DataTypes type;
      unsigned int components = 1;
      std::uint64_t offset = 0;
      Sections section = NO_SECTION;
    };

    // Type of global world coordinates
    using GlobalCoordinate = typename GC::GlobalCoordinate;

    // Template representing a grid-function that is created in getPointData() and getCellData()
    // with Context either Vtk::PointContext or Vek::CellContext, respectively.
    // To each GridCreator a GridFunction is associated, see, e.g. Vtk::ContinuousGridFunction
    // or Vtk::LagrangeGridFunction.
    template <class Context>
    using GridFunction = typename Vtk::AssociatedGridFunction<GC, FieldType, Context>::type;

  public:
    using GridCreator = GC;

    /// GridFunction representing the data stored on the points in the file
    using PointGridFunction = GridFunction<Vtk::PointContext>;

    /// GridFunction representing the data stored on the cells in the file
    using CellGridFunction = GridFunction<Vtk::CellContext>;

  public:
    /// Constructor. Creates a new GridCreator with the passed factory
    /**
     * \param args... Either pass a GridFactory by reference or shared_ptr, or a list of arguments
     *                passed to the constructor of a Dune::GridFactory (typically and empty parameter
     *                list). See the constructor of \ref GridCreatorInterface and the GridCreator
     *                passed to this reader.
     **/
    template <class... Args,
      std::enable_if_t<std::is_constructible<GridCreator, Args...>::value,int> = 0>
    explicit VtkReader (Args&&... args)
      : VtkReader(std::make_shared<GridCreator>(std::forward<Args>(args)...))
    {}

    /// Constructor. Stores the references in a non-destroying shared_ptr
    explicit VtkReader (GridCreator& creator)
      : VtkReader(stackobject_to_shared_ptr(creator))
    {}

    /// Constructor. Stores the shared_ptr
    explicit VtkReader (std::shared_ptr<GridCreator> creator)
      : creator_(std::move(creator))
    {}

    /// Read the grid from file with `filename` into the GridCreator
    /**
     * This function fills internal data containers representing the information from the
     * passed file.
     *
     * \param filename     The name of the input file
     * \param fillCreator  If `false`, only fill internal data structures, if `true`, pass
     *                     the internal data to the GridCreator. [true]
     **/
    void read (std::string const& filename, bool fillCreator = true);

    /// Obtains the creator of the reader
    GridCreator& gridCreator ()
    {
      return *creator_;
    }

    /// Construct the actual grid using the GridCreator
    /// [[expects: read_ == true]]
    std::unique_ptr<Grid> createGrid () const
    {
      return creator_->createGrid();
    }

    /// Construct a grid-function representing the point-data with the given name
    /// [[expects: read_ == true]]
    GridFunction<Vtk::PointContext> getPointData (std::string const& name) const
    {
      auto data_it = dataArray_.find("PointData." + name);
      auto point_it = pointData_.find("PointData." + name);
      VTK_ASSERT_MSG(data_it != dataArray_.end() && point_it != pointData_.end(),
        "The data to extract is not found in point-data. Try `getCellData()` instead!");
      VTK_ASSERT(data_it->second.section == POINT_DATA);

      return {*creator_, point_it->second, data_it->second.name, data_it->second.components,
              data_it->second.type, vec_types, vec_offsets, vec_connectivity};
    }

    /// Return a vector of DataArrayAttributes for all POINT_DATA blocks
    /// [[expects: read_ == true]]
    std::vector<DataArrayAttributes> getPointDataAttributes () const
    {
      std::vector<DataArrayAttributes> attributes;
      attributes.reserve(pointData_.size());
      for (auto const& da : dataArray_) {
        if (da.second.section == POINT_DATA)
          attributes.push_back(da.second);
      }
      return attributes;
    }

    /// Construct a grid-function representing the cell-data with the given name
    /// [[expects: read_ == true]]
    GridFunction<Vtk::CellContext> getCellData (std::string const& name) const
    {
      auto data_it = dataArray_.find("CellData." + name);
      auto cell_it = cellData_.find("CellData." + name);
      VTK_ASSERT_MSG(data_it != dataArray_.end() && cell_it != cellData_.end(),
        "The data to extract is not found in cell-data. Try `getPointData()` instead!");
      VTK_ASSERT(data_it->second.section == CELL_DATA);

      return {*creator_, cell_it->second, data_it->second.name, data_it->second.components,
              data_it->second.type, vec_types, vec_offsets, vec_connectivity};
    }

    /// Return a vector of DataArrayAttributes for all CELL_DATA blocks
    /// [[expects: read_ == true]]
    std::vector<DataArrayAttributes> getCellDataAttributes () const
    {
      std::vector<DataArrayAttributes> attributes;
      attributes.reserve(cellData_.size());
      for (auto const& da : dataArray_) {
        if (da.second.section == CELL_DATA)
          attributes.push_back(da.second);
      }
      return attributes;
    }

    /// Advanced read methods
    /// @{

    /// Read the grid from an input stream, referring to a .vtu file, into the GridFactory \ref factory_
    /**
     * \param input   A STL input stream to read the VTK file from.
     * \param create  If `false`, only fill internal data structures, if `true`, also create the grid. [true]
     **/
    void readSerialFileFromStream (std::ifstream& input, bool create = true);

    /// Read the grid from and input stream, referring to a .pvtu file, into the GridFactory \ref factory_
    /**
     * \param input   A STL input stream to read the VTK file from.
     * \param create  If `false`, only fill internal data structures, if `true`, also create the grid. [true]
     **/
    void readParallelFileFromStream (std::ifstream& input, int rank, int size, bool create = true);

    /// Insert all internal data to the GridCreator
    /// NOTE: requires an aforegoing call to \ref read()
    void fillGridCreator (bool insertPieces = true);

    /// @}

    /// Return the filenames of parallel pieces
    std::vector<std::string> const& pieces () const
    {
      return pieces_;
    }

#ifndef DOXYGEN
    // Implementation of the FileReader interface
    static void fillFactoryImpl (GridFactory<Grid>& factory, std::string const& filename)
    {
      VtkReader reader{factory};
      reader.read(filename);
    }
#endif

  private:
    // Read values stored on the cells with ID `id`
    Sections readCellData (std::ifstream& input, std::string id);

    template <class F, class H>
    void readCellDataAppended (MetaType<F>, MetaType<H>, std::ifstream& input, std::string id);

    // Read values stored on the points with ID `id`
    Sections readPointData (std::ifstream& input, std::string id);

    template <class F, class H>
    void readPointDataAppended (MetaType<F>, MetaType<H>, std::ifstream& input, std::string id);


    // Read vertex coordinates from `input` stream and store in into `factory`
    Sections readPoints (std::ifstream& input, std::string id);

    template <class F, class H>
    void readPointsAppended (MetaType<F>, MetaType<H>, std::ifstream& input);


    // Read cell type, cell offsets and connectivity from `input` stream
    Sections readCells (std::ifstream& input, std::string id);

    template <class H>
    void readCellsAppended (MetaType<H>, std::ifstream& input);

    // Read data from appended section in vtk file, starting from `offset`
    template <class FloatType, class HeaderType>
    void readAppended (std::ifstream& input, std::vector<FloatType>& values, HeaderType offset);

    // Test whether line belongs to section
    bool isSection (std::string line,
                    std::string key,
                    Sections current,
                    Sections parent = NO_SECTION) const
    {
      bool result = line.substr(1, key.length()) == key;
      if (result && current != parent)
        DUNE_THROW(Exception , "<" << key << "> in wrong section." );
      return result;
    }

    // Convert a section into a string
    std::string toString (Sections s) const;

    // Find beginning of appended binary data
    std::uint64_t findAppendedDataPosition (std::ifstream& input) const;

    // Read attributes from current xml tag
    std::map<std::string, std::string> parseXml (std::string const& line, bool& closed);

    // clear all vectors
    void clear ();

    auto comm () const
    {
      return MPIHelper::getCollectiveCommunication();
    }

  private:
    std::shared_ptr<GridCreator> creator_;

    /// Data format, i.e. ASCII, BINARY or COMPRESSED. Read from xml attributes.
    Vtk::FormatTypes format_;

    /// Type of compression algorithm used for binary data
    Vtk::CompressorTypes compressor_;

    // Temporary data to construct the grid elements
    std::vector<GlobalCoordinate> vec_points;
    std::vector<std::uint64_t> vec_point_ids;   //< Global unique vertex ID
    std::vector<std::uint8_t> vec_types;        //< VTK cell type ID
    std::vector<std::int64_t> vec_offsets;      //< offset of vertices of cell
    std::vector<std::int64_t> vec_connectivity; //< vertex indices of cell

    std::size_t numberOfCells_ = 0;   //< Number of cells in the grid
    std::size_t numberOfPoints_ = 0;  //< Number of vertices in the grid

    // offset information for appended data
    // map: id -> {Name,DataType,NumberOfComponents,Offset}
    std::map<std::string, DataArrayAttributes> dataArray_;

    // storage for internal point and cell data
    // map: id -> vector of data entries per point/cell
    std::map<std::string, std::vector<FieldType>> pointData_;
    std::map<std::string, std::vector<FieldType>> cellData_;

    // vector of filenames of parallel pieces
    std::vector<std::string> pieces_;

    /// Offset of beginning of appended data
    std::uint64_t offset0_ = 0;

    bool read_ = false;
  };

  // deduction guides
  template <class Grid>
  VtkReader (GridFactory<Grid>&)
    -> VtkReader<Grid, Vtk::ContinuousGridCreator<Grid>>;

  template <class GridCreator,
    class = std::void_t<typename GridCreator::Grid>>
  VtkReader (GridCreator&)
    -> VtkReader<typename GridCreator::Grid, GridCreator>;

  template <class GridCreator,
    class = std::void_t<typename GridCreator::Grid>>
  VtkReader (std::shared_ptr<GridCreator>)
    -> VtkReader<typename GridCreator::Grid, GridCreator>;

} // end namespace Dune

#include "vtkreader.impl.hh"
