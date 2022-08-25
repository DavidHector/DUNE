set(GRID_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../doc/")

dune_add_test(SOURCES vtkreader.cc
  LINK_LIBRARIES dunevtk)

dune_add_test(SOURCES legacyvtkwriter.cc
  LINK_LIBRARIES dunevtk)


dune_add_test(SOURCES vtkwriter.cc
  LINK_LIBRARIES dunevtk
  CMAKE_GUARD dune-functions_FOUND)

dune_add_test(SOURCES benchmark.cc
  LINK_LIBRARIES dunevtk
  CMAKE_GUARD dune-functions_FOUND)

dune_add_test(SOURCES datacollector.cc
  LINK_LIBRARIES dunevtk
  CMAKE_GUARD dune-functions_FOUND)

dune_add_test(SOURCES datareader.cc
  LINK_LIBRARIES dunevtk
  CMAKE_GUARD dune-functions_FOUND)

dune_add_test(SOURCES dunegridvtkreader.cc
  LINK_LIBRARIES dunevtk
  CMAKE_GUARD "dune-uggrid_FOUND OR dune-alugrid_FOUND")

dune_add_test(SOURCES lagrangepoints.cc
  LINK_LIBRARIES dunevtk)

dune_add_test(SOURCES lagrangereader.cc
  COMPILE_DEFINITIONS "GRID_PATH=\"${GRID_PATH}\""
  LINK_LIBRARIES dunevtk
  CMAKE_GUARD dune-foamgrid_FOUND)

dune_add_test(SOURCES structuredgridwriter.cc
  LINK_LIBRARIES dunevtk
  CMAKE_GUARD dune-functions_FOUND)

dune_add_test(SOURCES geometrygrid.cc
  LINK_LIBRARIES dunevtk
  CMAKE_GUARD dune-functions_FOUND)

dune_add_test(SOURCES timeserieswriter.cc
  LINK_LIBRARIES dunevtk
  CMAKE_GUARD dune-functions_FOUND)

dune_add_test(SOURCES pvdwriter.cc
  LINK_LIBRARIES dunevtk
  CMAKE_GUARD dune-functions_FOUND)

dune_add_test(SOURCES vectorwriter.cc
  LINK_LIBRARIES dunevtk
  CMAKE_GUARD dune-functions_FOUND)

if (dune-polygongrid_FOUND)
  # CMAKE_GUARD can not be used, since a dummy target is created and linked against dunepolygongrid
  dune_add_test(SOURCES polygongrid.cc
    LINK_LIBRARIES dunevtk dunepolygongrid)
endif ()

add_subdirectory(test)