// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PYTHON_VTK_WRITER_HH
#define DUNE_PYTHON_VTK_WRITER_HH

#include <dune/vtk/vtkwriter.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

namespace Dune
{

  namespace Vtk
  {

    // registerVTKWriter
    // -----------------

    template< class Writer, class... options >
    inline static void registerVtkWriter ( pybind11::handle scope,
                                           pybind11::class_< Writer, options... > cls )
    {
      using GridView = typename Writer::GridView;
      using VirtualizedGF = Dune::Vtk::Function<GridView>;

      cls.def( pybind11::init( [] ( GridView &grid,
               Vtk::FormatTypes format,
               Vtk::DataTypes datatype,
               Vtk::DataTypes headertype,
               pybind11::kwargs kwargs) {
        return new Writer( grid, format, datatype, headertype ); 
        }), pybind11::arg("grid"),
            pybind11::arg("format") = Vtk::FormatTypes::BINARY,
            pybind11::arg("datatype") = Vtk::DataTypes::FLOAT32,
            pybind11::arg("headertype") = Vtk::DataTypes::UINT32,
            pybind11::keep_alive< 1, 2 >());

      cls.def( "write",
          [] ( Writer &writer, const std::string &name ) {
            writer.write( name );
          },
          pybind11::arg("name") );
      cls.def( "write",
          [] ( Writer &writer, const std::string &name, int number ) {
            std::stringstream s; s << name << std::setw(5) << std::setfill('0') << number;
            writer.write( s.str());
          },
          pybind11::arg("name"),
          pybind11::arg("number") );

      cls.def( "addPointData",
          [] ( Writer &writer, VirtualizedGF &f,
               RangeTypes range, DataTypes data
             ) {
            f.setRangeType(range);
            f.setDataType(data);
            writer.addPointData(f);
          },
          pybind11::keep_alive< 1, 2 >(),
          pybind11::arg("f"),
          pybind11::arg("range")=RangeTypes::AUTO,
          pybind11::arg("data")=DataTypes::FLOAT32
        );
      cls.def( "addPointData",
          [] ( Writer &writer, VirtualizedGF &f,
               std::string &name,
               RangeTypes range, DataTypes data
             ) {
            f.setName(name);
            f.setRangeType(range);
            f.setDataType(data);
            writer.addPointData(f);
          },
          pybind11::keep_alive< 1, 2 >(),
          pybind11::arg("f"), pybind11::arg("name"),
          pybind11::arg("range")=RangeTypes::AUTO,
          pybind11::arg("data")=DataTypes::FLOAT32
        );
      cls.def( "addPointData",
          [] ( Writer &writer, VirtualizedGF &f,
               std::string &name, 
               std::vector<int> &components,
               RangeTypes range, DataTypes data
             ) {
            f.setName(name);
            f.setRangeType(range);
            f.setDataType(data);
            f.setComponents(components);
            writer.addPointData(f);
          },
          pybind11::keep_alive< 1, 2 >(),
          pybind11::arg("f"), pybind11::arg("name"),
          pybind11::arg("components"),
          pybind11::arg("range")=RangeTypes::AUTO,
          pybind11::arg("data")=DataTypes::FLOAT32
        );
      cls.def( "addPointData",
          [] ( Writer &writer, VirtualizedGF &f,
               FieldInfo &info) {
            f.setFieldInfo(info);
            writer.addPointData(f);
          },
          pybind11::keep_alive< 1, 2 >(),
          pybind11::arg("f"), pybind11::arg("info") );
      cls.def( "addPointData",
          [] ( Writer &writer, VirtualizedGF &f,
               std::vector<int> &components, FieldInfo &info ) {
            f.setFieldInfo(info);
            f.setComponents(components);
            writer.addPointData(f);
          },
          pybind11::keep_alive< 1, 2 >(),
          pybind11::arg("f"), pybind11::arg("components"), pybind11::arg("info") );

      cls.def( "addCellData",
          [] ( Writer &writer, VirtualizedGF &f,
               RangeTypes range, DataTypes data
             ) {
            f.setRangeType(range);
            f.setDataType(data);
            writer.addCellData(f);
          },
          pybind11::keep_alive< 1, 2 >(),
          pybind11::arg("f"),
          pybind11::arg("range")=RangeTypes::AUTO,
          pybind11::arg("data")=DataTypes::FLOAT32
        );
      cls.def( "addCellData",
          [] ( Writer &writer, VirtualizedGF &f,
               std::string &name,
               RangeTypes range, DataTypes data
             ) {
            f.setName(name);
            f.setRangeType(range);
            f.setDataType(data);
            writer.addCellData(f);
          },
          pybind11::keep_alive< 1, 2 >(),
          pybind11::arg("f"), pybind11::arg("name"),
          pybind11::arg("range")=RangeTypes::AUTO,
          pybind11::arg("data")=DataTypes::FLOAT32
        );
      cls.def( "addCellData",
          [] ( Writer &writer, VirtualizedGF &f,
               std::string &name, 
               std::vector<int> &components,
               RangeTypes range, DataTypes data
             ) {
            f.setName(name);
            f.setRangeType(range);
            f.setDataType(data);
            f.setComponents(components);
            writer.addCellData(f);
          },
          pybind11::keep_alive< 1, 2 >(),
          pybind11::arg("f"), pybind11::arg("name"),
          pybind11::arg("components"),
          pybind11::arg("range")=RangeTypes::AUTO,
          pybind11::arg("data")=DataTypes::FLOAT32
        );
      cls.def( "addCellData",
          [] ( Writer &writer, VirtualizedGF &f,
               FieldInfo &info) {
            f.setFieldInfo(info);
            writer.addCellData(f);
          },
          pybind11::keep_alive< 1, 2 >(),
          pybind11::arg("f"), pybind11::arg("info") );
      cls.def( "addCellData",
          [] ( Writer &writer, VirtualizedGF &f,
               std::vector<int> &components, FieldInfo &info ) {
            f.setFieldInfo(info);
            f.setComponents(components);
            writer.addCellData(f);
          },
          pybind11::keep_alive< 1, 2 >(),
          pybind11::arg("f"), pybind11::arg("components"), pybind11::arg("info") );
    }

  } // namespace Vtk

} // namespace Dune

#endif // DUNE_PYTHON_VTK_WRITER_HH
