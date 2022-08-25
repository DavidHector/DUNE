#install headers
install(FILES circularbuffer.hh conjugationfactor.hh linesearchcriterion.hh
  linesearch.hh linesearchproposal.hh nonlinopt.hh stoppingcriterion.hh
  unconstrained.hh userinterface.hh vectorclass.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/nonlinopt)
