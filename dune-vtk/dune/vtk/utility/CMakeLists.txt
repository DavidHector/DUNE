dune_add_library("filesystem" OBJECT
  filesystem.cc)

#install headers
install(FILES
  arguments.hh
  concepts.hh
  enum.hh
  errors.hh
  filesystem.hh
  lagrangepoints.hh
  lagrangepoints.impl.hh
  string.hh
  uid.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vtk/utility)

  add_subdirectory(test)