#pragma once

#include <memory>
#include <string>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridfactory.hh>

namespace Dune
{
  namespace Vtk
  {
    template <class Grid, class FilerReaderImp>
    class FileReader
    {
    private:
      // type of underlying implementation, for internal use only
      using Implementation = FilerReaderImp;

      /// \brief An accessor class to call protected members of reader implementations.
      struct Accessor : public Implementation
      {
        template <class... Args>
        static std::unique_ptr<Grid> createGridFromFileImpl (Args&&... args)
        {
          return Implementation::createGridFromFileImpl(std::forward<Args>(args)...);
        }

        template <class... Args>
        static void fillFactoryImpl (Args&&... args)
        {
          return Implementation::fillFactoryImpl(std::forward<Args>(args)...);
        }
      };

    public:
      /// Reads the grid from a file with filename and returns a unique_ptr to the created grid.
      /// Redirects to concrete implementation of derivated class.
      template <class... Args>
      static std::unique_ptr<Grid> createGridFromFile (const std::string &filename, Args&&... args)
      {
        return Accessor::createGridFromFileImpl(filename, std::forward<Args>(args)...);
      }

      /// Reads the grid from a file with filename into a grid-factory.
      /// Redirects to concrete implementation of derivated class.
      template <class... Args>
      static void fillFactory (GridFactory<Grid> &factory,
                              const std::string &filename,
                              Args&&... args)
      {
        Accessor::fillFactoryImpl(factory, filename, std::forward<Args>(args)...);
      }

    protected: // default implementations

      // Default implementation, redirects to fillFactory implementation.
      template <class... Args>
      static std::unique_ptr<Grid> createGridFromFileImpl (const std::string &filename,
                                                          Args&&... args)
      {
        GridFactory<Grid> factory;
        fillFactory(factory, filename, std::forward<Args>(args)...);

        return std::unique_ptr<Grid>{ factory.createGrid() };
      }

      // Default implementation for reading into grid-factory: produces a runtime-error.
      template <class... Args>
      static void fillFactoryImpl (GridFactory<Grid> &/*factory*/,
                                  const std::string &/*filename*/,
                                  Args&&... /*args*/)
      {
        DUNE_THROW(NotImplemented,
          "GridReader using a factory argument not implemented for concrete reader implementation.");
      }
    };

  } // end namespace Vtk
} // end namespace Dune
