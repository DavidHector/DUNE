#install headers
install(FILES
  vtkimagedatawriter.hh
  vtkimagedatawriter.impl.hh
  vtkrectilineargridwriter.hh
  vtkrectilineargridwriter.impl.hh
  vtkstructuredgridwriter.hh
  vtkstructuredgridwriter.impl.hh
  vtkunstructuredgridwriter.hh
  vtkunstructuredgridwriter.impl.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/vtk/writers)
