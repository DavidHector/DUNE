cmake_minimum_required(VERSION 2.8.12)
project(dune-nonlinopt CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common)

# configuration as a dune module
if(dune-common_FOUND)
  message("dune-common found, configuring as module")

  list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
    ${dune-common_MODULE_PATH})

  #include the dune macros
  include(DuneMacros)

  find_package(Eigen3)
  if (EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
    add_compile_definitions(PUBLIC HAVE_EIGEN)
  endif()
  find_package(muParser)
  if (muParser_FOUND)
    link_libraries(${MUPARSER_LIBRARY})
    include_directories(PUBLIC ${MUPARSER_INCLUDE_DIR})
    add_compile_definitions(PUBLIC HAVE_MUPARSER)
  endif()


  # start a dune project with information from dune.module
  dune_project()
  dune_enable_all_packages()

  add_subdirectory("dune")
  add_subdirectory("examples")
  add_subdirectory("cmake/modules")
  add_subdirectory("doc")
  add_subdirectory("test")

  # finalize the dune project, e.g. generating config.h etc.
  finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

#configuration as a standalone package
else()
  message("dune-common not found, configuring standalone package")
  # specify the C++ standard
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED True)

  # provide option to specify non-standard software installs
  set(CMAKE_PREFIX_PATH "" CACHE
      STRING "Path to software in non-standard locations.")

  # search for external dependencies
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)
  find_package(Eigen3)
  if (EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
    add_compile_definitions(PUBLIC HAVE_EIGEN)
  endif()
  find_package(muParser)
  if (muParser_FOUND)
    link_libraries(${MUPARSER_LIBRARY})
    include_directories(PUBLIC ${MUPARSER_INCLUDE_DIR})
    add_compile_definitions(PUBLIC HAVE_MUPARSER)
  endif()
  find_package(Doxygen)
  
  include_directories(${CMAKE_SOURCE_DIR})

  add_subdirectory("dune")
  add_subdirectory("examples")
  add_subdirectory("cmake/modules")

  # add documentation directory
  option(BUILD_DOCUMENTATION "Create Doxygen documentation." ON)
  if (BUILD_DOCUMENTATION AND Doxygen_FOUND)
    add_subdirectory("doc")
  endif()

  if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)

    #add test directory
    if (BUILD_TESTING)
      add_subdirectory("test")
    endif()
  endif()
endif()
