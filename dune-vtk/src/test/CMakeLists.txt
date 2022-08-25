dune_add_test(SOURCES reader_writer_test.cc
              LINK_LIBRARIES dunevtk)

dune_add_test(SOURCES parallel_reader_writer_test.cc
              LINK_LIBRARIES dunevtk)

dune_add_test(SOURCES mixed_element_test.cc
              LINK_LIBRARIES dunevtk)
