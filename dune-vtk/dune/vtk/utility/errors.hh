#pragma once

#include <dune/common/exceptions.hh>

/**
 * \file
 * \brief Macro for wrapping error checks and throwing exceptions
 */

namespace Dune {

class VtkError : public Exception {};

}

/**
 * \brief check if condition \a cond holds; otherwise, throw a VtkError with a message.
 */
#define VTK_ASSERT_MSG(cond, text)      \
  do {                                  \
    if (!(cond))                        \
      DUNE_THROW(Dune::VtkError, text); \
  } while (false)


/**
 * \brief check if condition \a cond holds; otherwise, throw a VtkError.
 */
#define VTK_ASSERT(cond)                \
  do {                                  \
    if (!(cond))                        \
      DUNE_THROW(Dune::VtkError, #cond); \
  } while (false)
