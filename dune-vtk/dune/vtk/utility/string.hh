#pragma once

#include <algorithm>
#include <cctype>
#include <locale>
#include <sstream>
#include <string>

namespace Dune
{
  namespace Vtk
  {
    /// convert all characters in a string to upper case
    inline std::string to_upper(std::string input)
    {
      for (auto& c : input)
        c = toupper(c);
      return input;
    }

    /// convert all characters in a string to upper case
    inline std::string to_lower(std::string input)
    {
      for (auto& c : input)
        c = tolower(c);
      return input;
    }

    /// trim a string from the left
    inline std::string& ltrim(std::string& str)
    {
      auto it =  std::find_if(str.begin(), str.end(), [](char ch)
      {
        return !std::isspace<char>(ch, std::locale::classic());
      });
      str.erase(str.begin() , it);
      return str;
    }

    /// trim a string from the right
    inline std::string& rtrim(std::string& str)
    {
      auto it =  std::find_if(str.rbegin(), str.rend(), [](char ch)
      {
        return !std::isspace<char>(ch, std::locale::classic());
      });
      str.erase(it.base(), str.end());
      return str;
    }

    /// trim a string from both sides
    inline std::string& trim(std::string& str)
    {
      return ltrim(rtrim(str));
    }

    /// trim a (copy of the) string from both sides
    inline std::string trim_copy(std::string const& str)
    {
      auto s = str;
      return trim(s);
    }


    template <class InputIter, class T, class Func>
    void split(InputIter first, InputIter end, T const& t, Func f)
    {
      if (first == end)
        return;

      while (true) {
        InputIter found = std::find(first, end, t);
        f(first, found);
        if (found == end)
          break;
        first = ++found;
      }
    }

    template <class InputIter, class SeparatorIter, class Func>
    void split(InputIter first, InputIter end, SeparatorIter s_first, SeparatorIter s_end, Func f)
    {
      if (first == end)
        return;

      while (true) {
        InputIter found = std::find_first_of(first, end, s_first, s_end);
        f(first, found);
        if (found == end)
          break;
        first = ++found;
      }
    }

    /// Replace all occurences of substring `from` with `to` in source `str`.
    inline void replaceAll(std::string& str, std::string const& from, std::string const& to)
    {
      if (from.empty())
        return;
      std::size_t start_pos = 0;
      while ((start_pos = str.find(from, start_pos)) != std::string::npos)
      {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
      }
    }


    template <class InputIter>
    std::string join (InputIter first, InputIter end, std::string sep = " ")
    {
      if (first == end)
        return "";

      std::ostringstream os;
      os << *first++;
      while (first != end)
        os << sep << *first++;
      return os.str();
    }

  } // end namespace Vtk
} // end namspace Dune
