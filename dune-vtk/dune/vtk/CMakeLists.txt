dune_add_library("vtktypes" OBJECT
  types.cc)

#install headers
install(FILES
  datacollectorinterface.hh
  datacollectorinterface.impl.hh
  defaultvtkfunction.hh
  filereader.hh
  filewriter.hh
  gridcreatorinterface.hh
  legacyvtkfunction.hh
  pvdwriter.hh
  pvdwriter.impl.hh
  function.hh
  localfunction.hh
  localfunctioninterface.hh
  vtkreader.hh
  vtkreader.impl.hh
  vtktimeserieswriter.hh
  vtktimeserieswriter.impl.hh
  types.hh
  vtkwriter.hh
  vtkwriterinterface.hh
  vtkwriterinterface.impl.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vtk)

add_subdirectory(datacollectors)
add_subdirectory(gridcreators)
add_subdirectory(gridfunctions)
add_subdirectory(test)
add_subdirectory(utility)
add_subdirectory(writers)
