#pragma once

#include <type_traits>

namespace Dune
{
  namespace Vtk
  {

    template <class E, class Integer,
      std::enable_if_t<std::is_enum<E>::value, int> = 0>
    constexpr bool is_a(E a, Integer b)
    {
      return (int(a) & int(b)) != 0;
    }

  } // end namespace Vtk
} // end namespace Dune
