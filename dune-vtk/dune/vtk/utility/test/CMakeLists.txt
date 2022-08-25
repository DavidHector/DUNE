dune_add_test(SOURCES test-arguments.cc)

dune_add_test(NAME test-lagrange1
              SOURCES test-lagrange.cc
              COMPILE_DEFINITIONS "CHECKDIM=1")

dune_add_test(NAME test-lagrange2
              SOURCES test-lagrange.cc
              COMPILE_DEFINITIONS "CHECKDIM=2")

dune_add_test(NAME test-lagrange3
              SOURCES test-lagrange.cc
              COMPILE_DEFINITIONS "CHECKDIM=3")
