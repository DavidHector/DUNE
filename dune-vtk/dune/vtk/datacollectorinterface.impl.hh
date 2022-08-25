#pragma once

#include <dune/geometry/referenceelements.hh>

namespace Dune {
namespace Vtk {

template <class GV, class D, class P>
  template <class T, class VtkFunction>
std::vector<T> DataCollectorInterface<GV,D,P>
  ::cellDataImpl (VtkFunction const& fct) const
{
  std::vector<T> data;
  data.reserve(this->numCells() * fct.numComponents());

  auto localFct = localFunction(fct);
  for (auto const& e : elements(gridView_, partition)) {
    localFct.bind(e);
    auto refElem = referenceElement<T,dim>(e.type());
    for (int comp = 0; comp < fct.numComponents(); ++comp)
      data.emplace_back(localFct.evaluate(comp, refElem.position(0,0)));
    localFct.unbind();
  }
  return data;
}

} } // end namespace Dune::Vtk
