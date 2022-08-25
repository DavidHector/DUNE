#pragma once

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>

namespace Dune
{
  namespace Vtk
  {
    inline std::string uid (std::size_t len = 8)
    {
      static const auto digits = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
      static const int N = std::strlen(digits);

      std::string id(len,' ');
      for (std::size_t i = 0; i < len; ++i)
        id[i] = digits[std::rand()%N];

      return id;
    }

  } // end namespace Vtk
} // end namespace Dune
