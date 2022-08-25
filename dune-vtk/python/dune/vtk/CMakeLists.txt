add_python_targets(vtk
  __init__
)
dune_add_pybind11_module(NAME _vtk)
set_property(TARGET _vtk PROPERTY LINK_LIBRARIES dunecommon dunegeometry dunegrid dunevtk APPEND)
