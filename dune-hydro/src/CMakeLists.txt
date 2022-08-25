add_executable("firstexample" firstexample.cc)
target_link_dune_default_libraries(firstexample)
target_link_libraries(firstexample ${GDAL_LIBRARIES})
target_include_directories(firstexample PUBLIC ${GDAL_INCLUDE_DIRS})

add_executable("secondexample" secondexample.cc)
target_link_dune_default_libraries(secondexample)
target_link_libraries(secondexample ${GDAL_LIBRARIES})
target_include_directories(secondexample PUBLIC ${GDAL_INCLUDE_DIRS})

add_executable("thirdexample" thirdexample.cc)
target_link_dune_default_libraries(thirdexample)
target_link_libraries(thirdexample ${GDAL_LIBRARIES})
target_include_directories(thirdexample PUBLIC ${GDAL_INCLUDE_DIRS})

add_executable("fourthexample" fourthexample.cc)
target_link_dune_default_libraries(fourthexample)
target_link_libraries(fourthexample  GDAL::GDAL)

add_executable("readtiff" readtiff.cc)
target_link_dune_default_libraries(readtiff)
target_link_libraries(readtiff GDAL::GDAL)
target_link_libraries(readtiff NetCDF::NetCDF)

add_executable("bw90mpre" bw90mpre.cc)
target_link_dune_default_libraries(bw90mpre)
target_link_libraries(bw90mpre GDAL::GDAL)
target_link_libraries(bw90mpre NetCDF::NetCDF)

add_executable("stok" stok.cc)
target_link_dune_default_libraries(stok)
target_link_libraries(stok GDAL::GDAL)
target_link_libraries(stok NetCDF::NetCDF)

add_executable("stok_long_term" stok_long_term.cc)
target_link_dune_default_libraries(stok_long_term)
target_link_libraries(stok_long_term GDAL::GDAL)
target_link_libraries(stok_long_term NetCDF::NetCDF)

add_executable("coupledexample" coupledexample.cc)
target_link_dune_default_libraries(coupledexample)
target_link_libraries(coupledexample GDAL::GDAL)
target_link_libraries(coupledexample NetCDF::NetCDF)

dune_symlink_to_source_files(FILES fourthexample.ini thirdexample.ini secondexample.ini firstexample.ini readtiff.ini bw90mpre.ini stok.ini coupledexample.ini stok_long_term.ini)

add_executable("dune-hydro" dune-hydro.cc)
target_link_dune_default_libraries(dune-hydro)
target_link_libraries(dune-hydro ${GDAL_LIBRARIES})
target_include_directories(dune-hydro PUBLIC ${GDAL_INCLUDE_DIRS})
