#pragma once

#include <string>
#include <vector>

#include "string.hh"

namespace Dune
{
  namespace Vtk
  {
    // A minimalistic filesystem class
    class Path
        : public std::vector<std::string>
    {
      using Super = std::vector<std::string>;
      using iterator = Super::iterator;
      using const_iterator = Super::const_iterator;

    public:
#ifdef _WIN32
      static constexpr char preferredSeparator = '\\';
#else
      static constexpr char preferredSeparator = '/';
#endif

    public:
      Path() = default;

      // NOTE: implicit conversion is allowed here
      template <class String>
      Path(String const& p)
        : original(p)
      {
        split(p);
      }

      template <class InputIt>
      Path(InputIt it, InputIt end_it)
        : Super(it, end_it)
      {
        original = this->string();
      }

      template <class String>
      Path(std::initializer_list<String> const& list)
        : Path(list.begin(), list.end())
      {}

      /// Removes filename path component
      Path& removeFilename()
      {
        this->pop_back();
        return *this;
      }

      /// Returns the path of the parent path
      Path parentPath() const
      {
        return empty() ? Path() : Path(begin(), --end());
      }

      /// Returns filename path component
      Path filename() const
      {
        return empty() ? Path() : Path(back());
      }

      /// Returns the stem path component
      Path stem() const;

      /// Returns the file extension path component
      Path extension() const;

      /// Return the path as string
      std::string string() const;

      /// \brief Return whether a path is an absolute path.
      /** In Linux, test whether the path starts with `/`, in Windows whether it starts
        * with `[a-z]:\\`.
        **/
      static bool isAbsolute(std::string p);

      bool isAbsolute() const { return isAbsolute(original); }

      bool isRelative() const { return !isAbsolute(); }

      /// Check whether path is a regular file
      bool isFile() const;

      /// Check whether path is a regular file
      bool isDirectory() const;

      /// Lexicographically compares two paths
      bool operator==(Path const& p)
      {
        return this->string() == p.string();
      }

      /// Appends elements to the path
      Path& operator/=(Path const& p);

      /// output of the path
      template <class CharT, class Traits>
      friend std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits>& out, Path const& p)
      {
        out << '"' << p.string() << '"';
        return out;
      }

    protected:

      // split the path string into names separated by a `/`, remove relative directories,
      // like `.` or `..`, if possible.
      void split(std::string p);

    private:
      std::string original = "";
    };

    /// Test whether the path is a valid (existing and accessible) file / directory
    bool exists(Path const&);

    /// Create directory and non existing parent directories.
    bool createDirectories(Path const&);

    /// Returns the current path
    Path currentPath();

    /// Find the path of `a` relative to directory of `b`
    Path relative(Path const& a, Path const& b);

  } // end namespace Vtk
} // end namespace Dune
