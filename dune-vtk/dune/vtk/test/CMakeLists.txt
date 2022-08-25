dune_add_test(SOURCES test-map-datatypes.cc
  LINK_LIBRARIES dunevtk)

dune_add_test(SOURCES test-function.cc
  LINK_LIBRARIES dunevtk)

dune_add_test(SOURCES test-typededuction.cc
  LINK_LIBRARIES dunevtk)

dune_add_test(SOURCES test-vtkwriter.cc
  LINK_LIBRARIES dunevtk)