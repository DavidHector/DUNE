#install headers
install(FILES
  continuousdatacollector.hh
  discontinuousdatacollector.hh
  lagrangedatacollector.hh
  quadraticdatacollector.hh
  spdatacollector.hh
  structureddatacollector.hh
  unstructureddatacollector.hh
  yaspdatacollector.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vtk/datacollectors)

add_subdirectory(test)