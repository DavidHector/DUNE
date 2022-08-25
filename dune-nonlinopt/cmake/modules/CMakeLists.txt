if(dune-common_FOUND)
  set(modules "DuneNonlinoptMacros.cmake")

  install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
endif()
