if(dune-common_FOUND)
  # shortcut for creating the Doxyfile.in and Doxyfile
  add_doxygen_target()
else()
  # Find all the public headers
  file(GLOB_RECURSE PUBLIC_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/../../dune/nonlinopt/*.hh)

  set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/dune/nonlinopt)
  set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})
  set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
  set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  #Replace variables inside @@ with the current values
  configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

  file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

  add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                     DEPENDS ${PUBLIC_HEADERS}
                     COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                     COMMENT "Generating docs")

  add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})
endif()
