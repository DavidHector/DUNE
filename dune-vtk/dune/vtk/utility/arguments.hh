#pragma once

#include <type_traits>
#include <utility>

namespace Dune
{
  namespace Vtk
  {
    /// Extract the first argument of the variadic list that is the same as one of
    /// the types {T,...} and return its value.
    /**
     * This utility can be used to implement function parameters with flexible order.
     *
     * Example:
     * ```
     * template <class... Args>
     * void foo(Args const&... args)
     * {
     *   int i = getArg<int>(args...);
     *   double d = getArg<double>(args..., 42.0); // with default value
     * }
     * ```
     *
     * The arguments are tested from first to last. Thus, a default parmaeter could be given at
     * end of the variadic list that is chosen if no argument matches the requested type.
     **/
    template <class... T, class Arg0, class... Args>
    decltype(auto) getArg(Arg0&& arg0, Args&&... args)
    {
      using A = std::decay_t<Arg0>;
      if constexpr ((std::is_same_v<A,T> ||...))
        return std::forward<Arg0>(arg0);
      else
        return getArg<T...>(std::forward<Args>(args)...);
    }

  } // end namespace Vtk
} // end namespace Dune
