# OBJECT is needed since CMake 3.0.0
# but it conflicts with shared libs
set(_OBJECT_FLAG "OBJECT")
if(BUILD_SHARED_LIBS)
  set(_OBJECT_FLAG "")
endif()

dune_add_library(dunevtk
  _DUNE_TARGET_OBJECTS:filesystem_
  _DUNE_TARGET_OBJECTS:vtktypes_
  ADD_LIBS ${DUNE_LIBS})
