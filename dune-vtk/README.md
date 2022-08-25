# Dune-Vtk
File reader and writer for the VTK Format

## Summary
Provides structured and unstructured file writers for the VTK XML File Formats
that can be opened in the popular ParaView visualization application. Additionally
a file reader is provided to import VTK files into Dune grid and data objects.

## Installation Instructions
`dune-vtk` requires the DUNE core modules, version 2.7 or later.
Please see the [general instructions for building DUNE modules](https://www.dune-project.org/doc/installation)
for detailed instructions on how to build the module.

## Usage
The VTK writer works similar to the dune-grid `VTKWriter`. It needs to be bound
to a GridView and then data can be added to the points or cells in the grid.
Points are not necessarily grid vertices, but any coordinates placed inside the
grid cells, so the data must be provided as GridViewFunction to allow the local
evaluation in arbitrary local coordinates.

General interface of a VtkWriter
```c++
template <class GridView, class DataCollector = DefaultDataCollector<GridView>>
class Vtk[Type]Writer
{
public:
  // Constructor
  Vtk[Type]Writer(GridView, Vtk::FormatTypes = Vtk::FormatTypes::BINARY, Vtk::DataTypes = Vtk::DataTypes::FLOAT32);

  // Bind data to the writer
  Vtk[Type]Writer& addPointData(Function [, std::string name, int numComponents, Vtk::FormatTypes]);
  Vtk[Type]Writer& addCellData(Function [, std::string name, int numComponents, Vtk::FormatTypes]);

  // Write file with filename
  void write(std::string filename);
};
```
where `Function` is either a `GridViewFunction`, i.e. supports `bind()`, `unbind()`, and `localFunction(Function)`, or is a legacy `VTKFunction` from Dune-Grid. The optional parameters `name`, `numComponents` and `format` may be given for a `GridViewFunction`.

The parameter `Vtk::FormatTypes` is one of `Vtk::FormatTypes::ASCII`, `Vtk::FormatTypes::BINARY`, or `Vtk::FormatTypes::COMPRESSED` and `Vtk::DataTypes` is one of `Vtk::DataTypes::FLOAT32`, or `Vtk::DataTypes::FLOAT64`. The `[Type]` of a VtkWriter is one of `UnstructuredGrid`, `StructuredGrid`, `RectilinearGrid`, `ImageData`, or `Timeseries`, see below for details. A `DataCollector` may be specified to control how point and cell values are extracted from the `GridView` and the bound data. See `dune/vtk/datacollectors/` of a list of poissible types. The default datacollector extracts a connected grid with continuous data, where points are grid vertices.

See also the `src/` directory for more examples.

## Comparison with Dune::VTKWriter
In Dune-Grid there is a VTK writer available, that is a bit different from the
proposed one. A comparison:

| **Property**       | **Dune-Grid** | **Dune-Vtk** |
| ------------------ | :-----------: | :----------: |
| VTK version        | 0.1           | 0.1/1.0      |
| UnstructuredGrid   | **x**         | **x**        |
| PolyData           | (1d)          | -            |
| StructuredGrid     | -             | **x**        |
| RectilinearGrid    | -             | **x**        |
| ImageData          | -             | **x**        |
| ASCII              | **x**         | **x**        |
| BASE64             | **x**         | -            |
| APPENDED_RAW       | **x**         | **x**        |
| APPENDED_BASE64    | **x**         | -            |
| BASE64_COMPRESSED  | -             | -            |
| APPENDED_COMPRESSED| -             | **x**        |
| Parallel files     | **x**         | **x**        |
| Conforming Data    | **x**         | **x**        |
| NonConforming Data | **x**         | **x**        |
| Quadratic Data     | -             | **x**        |
| Higher-Order Data  | -             | **x**        |
| Subdivided Data    | **x**         | -            |
| Sequence (PVD)     | **x**         | **x**        |
| Timeseries         | -             | **x**        |

## Writers and Readers
Dune-Vtk provides nearly all file formats specified in VTK + 2 time series formats:
PVD and VTK-Timeseries.

### VtkUnstructuredGridWriter
Implements a VTK file format for unstructured grids with arbitrary element types
in 1d, 2d, and 3d. Coordinates are specified explicitly and a connectivity table +
element types are specified for all grid elements (of codim 0). Can be used with
all Dune grid types.

### VtkStructuredGridWriter
Implements a writer for grid composed of cube elements (lines, pixels, voxels) with
local numbering similar to Dunes `cube(d)` numbering. The coordinates of the vertices
can be arbitrary but the connectivity is implicitly given and equals that of
`Dune::YaspGrid` or `Dune::SPGrid`. Might be chosen as writer for a transformed
structured grid, using, e.g., a `GeometryGrid` meta-grid. See `src/geometrygrid.cc`
for an example.

### VtkRectilinearGridWriter
Rectilinear grids are tensor-product grids with given coordinates along the x, y,
and z axes. Therefore, the grid must allow to extract these 1d coordinates and a
specialization for a `StructuredDataCollector` must be provided, that implements
the `ordinates()` function. By default, it assumes constant grid spacing starting
from a lower left corner. For `YaspGrid` a specialization is implemented if the
coordinates type is `TensorProductCoordinates`. See `src/structuredgridwriter.cc`
for an example.

### VtkImageDataWriter
The *most structured* grid is composed of axis-parallel cube elements with constant
size along each axis. The is implemented in the VtkImageDataWriter. A specialization
of the `StructuredDataCollector` must implement `origin()` for the lower left corner,
`wholeExtent()` for the range of cell numbers along each axis in the global grid,
`extent()` for the range in the local grid, and `spacing()` for the constant grid
spacing in each direction.

### PvdWriter
A sequence writer, i.e. a collection of timestep files, in the ParaView Data (PVD)
format. Supports all VtkWriters for the timestep output. In each timestep a collection
(.pvd) file is created.

### VtkTimseriesWriter
A timeseries is a collection of timesteps stored in one file, instead of separate
files for each timestep value. Since in the `Vtk::FormatTypes::APPENDED` mode, the data is written
as binary blocks in the appended section of the file and references by an offset
in the XML DataArray attributes, it allows to reuse written data. An example of
usage is when the grid points and cells do not change over time, but just the
point-/cell-data. Then, the grid is written only once and the data is just appended.

Timeseries file are create a bit differently from other Vtk file. There, in the
first write the grid points and cells are stored in a separate file, and in each
timestep just the data is written also to temporary files. When you need the timeseries
file, these stored temporaries are collected and combined to one VTK file. Thus,
only the minimum amount of data is written in each timestep. The intermediate files
are stored, by default, in a `/tmp` folder, with (hopefully) fast write access.

## VtkReader
Reading unstructured grid files (.vtu files) and creating a new grid, using a GridFactory,
can be performed using the `VtkReader` class. The reader allows to create the grid in
multiple ways, by providing a `GridCreator` template parameter. The `ContinuousGridCreator`
reads the connectivity of the grid as it is and assumes that the elements are already
connected correctly. On the other hand, a `DiscontinuousGridCreator` reconnects separated
elements, by identifying matching coordinates of the cell vertices. See more possible
grid-creators in the directory `dune/vtk/gridcreators`.

General interface of a VtkReader
```c++
template <class Grid, class GridCreator = ContinuousGridCreator<Grid>, class FieldType = double>
class VtkReader
{
public:
  // Constructors
  VtkReader();                    // Construct a GridCreator with internal stored GridFactory
  VtkReader(GridFactory<Grid>&);  // Construct a GridCreator referencing the passed GridFactory
  VtkReader(GridCreator&);        // Reference the passed GridCreator

  // Read the data from a file with filename
  void read(std::string filename);

  // Construct the Grid from the data read before.
  std::unique_ptr<Grid> createGrid() const;

  // Static method to construct the Grid directly
  static std::unique_ptr<Grid> createGridFromFile(std::string file);

  // Extract data from the reader
  _PointDataGridFunction_ getPointData(std::string name) const;
  _CellDataGridFunction_  getCellData(std::string name) const;
};
```
where `Grid` is a dune grid type providing a `GridFactory` specialization, `GridCreator is
the policy type implementing how the raw data from the file is transformed in data that can
be passed to the GridFactory and `FieldType` is the data for for internal storage of value
data from the file, i.e. point-data or cell-data.

The grid can either be created using the static method `createGridFromFile()` or by first
constructing the `VtkReader` with its `GridCreator`, calling `read()` and finally `createGrid()`.
The latter allows to access data stored in the reader, like point-data or cell-data.

Value from point-data or cell-data cannot be accessed directly, but through the interface
of a grid-function. These grid-functions are provided by the `getPointData()` or `getCellData()`
member functions of the reader. The interface of a dune-functions grid-function concept is
implemented by these two types. The reason why the reader does not provide the data directly
is, that it is quiet complicated to associate the specific value to a DOF in the grid, since
the GridFactory is allows to change the global indexing and even to change to local indexing
in the elements such that even the local element coordinates might need a transformation
compared to that of the element stored in the file. All these renumbering and coordinate
transformations are performed by the grid-functions internally.

The VtkReader supports grid creation in parallel. If a partition file .pvtu is
provided, all partitions can be read by either one processor and distributed later on
(`SerialGridCreator`) or read directly in parallel (`ParallelGridCreator`). The later
is currently only available in dune-alugrid 2.6.
