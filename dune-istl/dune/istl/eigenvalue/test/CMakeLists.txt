add_executable(poweriterationtest cond2test.cc)
add_executable(mytest evtest.cc)
target_link_libraries(poweriterationtest PRIVATE dunecommon)
target_link_libraries(mytest PRIVATE dunecommon)
dune_add_test(
  TARGET poweriterationtest
  CMD_ARGS 40)

if(SuperLU_FOUND)
  add_executable(poweriterationsuperlutest cond2test.cc)
  add_executable(mysuperlutest evtest.cc)
  add_dune_superlu_flags(poweriterationsuperlutest)
  add_dune_superlu_flags(mysuperlutest)
  dune_add_test(
    TARGET poweriterationsuperlutest
    CMD_ARGS 40)
  target_link_libraries(poweriterationsuperlutest PRIVATE dunecommon)
  target_link_libraries(mysuperlutest PRIVATE dunecommon)
endif()

if(ARPACKPP_FOUND)
  add_executable(arpackpptest cond2test.cc)
  add_dune_arpackpp_flags(arpackpptest)
  target_link_libraries(arpackpptest PRIVATE dunecommon)
  dune_add_test(TARGET arpackpptest)
endif()

if(ARPACKPP_FOUND AND SuperLU_FOUND)
  dune_add_test(NAME arpackppsuperlutest
                SOURCES cond2test.cc)
endif()


