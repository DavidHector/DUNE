#install headers
install(FILES
  common.hh
  continuousgridfunction.hh
  lagrangegridfunction.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vtk/gridfunctions)

add_subdirectory(test)