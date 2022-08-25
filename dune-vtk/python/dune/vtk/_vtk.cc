// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <dune/vtk/types.hh>

#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( _vtk, module )
{
  {
    using namespace Dune::Vtk;
    pybind11::enum_< FormatTypes > formatTypes( module, "FormatTypes" );
    for (const auto e : formatTypesList)
      formatTypes.value( to_string(e).c_str(), e );
    pybind11::enum_< RangeTypes > rangeTypes( module, "RangeTypes" );
    for (const auto e : rangeTypesList)
      rangeTypes.value( to_string(e).c_str(), e );
    pybind11::enum_< DataTypes > dataTypes( module, "DataTypes" );
    for (const auto e : dataTypesLists)
      dataTypes.value( to_string(e).c_str(), e );

    pybind11::class_<FieldInfo> fieldInfo(module, "FieldInfo" );
    fieldInfo.def( pybind11::init( 
      [](std::string name, RangeTypes range, DataTypes data)
      { return new FieldInfo(name, range, data); }
      ), pybind11::arg("name"),
         pybind11::arg("range")=RangeTypes::AUTO,
         pybind11::arg("data")=DataTypes::FLOAT32 );
  }
}
