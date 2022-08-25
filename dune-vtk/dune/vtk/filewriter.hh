#pragma once

#include <optional>
#include <string>

namespace Dune
{
  namespace Vtk
  {
    class FileWriter
    {
    public:
      /// Virtual destructor
      virtual ~FileWriter () = default;

      /// Write to file given by `filename` and (optionally) store additional data in `dataDir`
      virtual std::string write (std::string const& filename, std::optional<std::string> dataDir = {}) const = 0;
    };

  } // end namespace Vtk
} // end namespace Dune
