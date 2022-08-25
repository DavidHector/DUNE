#pragma once

#include <numeric>
#include <type_traits>

#include <dune/common/typetraits.hh>
#include <dune/common/version.hh>

#include <dune/vtk/localfunction.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/utility/arguments.hh>
#include <dune/vtk/utility/concepts.hh>

namespace Dune
{
  // forward declarations
  template <class T, int N>
  class FieldVector;

  template <class T, int N, int M>
  class FieldMatrix;

  namespace Vtk
  {
    /// Wrapper class for functions allowing local evaluations.
    template <class GridView>
    class Function
    {
      using Element = typename GridView::template Codim<0>::Entity;
      using LocalDomain = typename Element::Geometry::LocalCoordinate;

      template <class F, class D>
      using Range = std::decay_t<std::result_of_t<F(D)>>;

    private:

      template <class T, int N>
      static auto sizeOfImpl (FieldVector<T,N>) -> std::integral_constant<int, N> { return {}; }

      template <class T, int N, int M>
      static auto sizeOfImpl (FieldMatrix<T,N,M>) -> std::integral_constant<int, N*M> { return {}; }

      static auto sizeOfImpl (...) -> std::integral_constant<int, 1> { return {}; }

      template <class T>
      static constexpr int sizeOf () { return decltype(sizeOfImpl(std::declval<T>()))::value; }

      static std::vector<int> allComponents(int n)
      {
        std::vector<int> components(n);
        std::iota(components.begin(), components.end(), 0);
        return components;
      }

    public:
      /// (1) Construct from a LocalFunction directly
      /**
      * \param localFct    A local-function, providing a `bind(Element)` and an `operator()(LocalDomain)`
      * \param name        The name to use as identification in the VTK file
      * \param components  A vector of component indices to extract from the range type
      * \param category    The \ref Vtk::RangeTypes category for the range. [Vtk::RangeTypes::UNSPECIFIED]
      * \param dataType    The \ref Vtk::DataTypes used in the output. [Vtk::DataTypes::FLOAT64]
      *
      * The arguments `category` and `dataType` can be passed in any order.
      *
      * NOTE: Stores the localFunction by value.
      **/
      template <class LF, class... Args,
        IsLocalFunction<LF,Element> = true>
      Function (LF&& localFct, std::string name, std::vector<int> components, Args&&... args)
        : localFct_(std::forward<LF>(localFct))
        , name_(std::move(name))
      {
        setComponents(std::move(components));
        setRangeType(getArg<Vtk::RangeTypes>(args..., Vtk::RangeTypes::UNSPECIFIED), components_.size());
        setDataType(getArg<Vtk::DataTypes>(args..., Vtk::DataTypes::FLOAT64));
      }

      /// (2) Construct from a LocalFunction directly
      /**
      * \param localFct   A local-function, providing a `bind(Element)` and an `operator()(LocalDomain)`
      * \param name    The name to use as identification in the VTK file
      * \param ncomps  Number of components of the pointwise data. Is extracted
      *                from the range type of the GridFunction if not given.
      *
      * Forwards all the other parmeters to the constructor (1)
      *
      * NOTE: Stores the localFunction by value.
      **/
      template <class LF, class... Args,
        IsLocalFunction<LF,Element> = true>
      Function (LF&& localFct, std::string name, int ncomps, Args&&... args)
        : Function(std::forward<LF>(localFct), std::move(name), allComponents(ncomps),
                   std::forward<Args>(args)...)
      {}

      /// (3) Construct from a LocalFunction directly.
      /**
       * Same as Constructor (1) or (2) but deduces the number of components from
       * the static range type of the local-function. This defaults to 1 of no
       * static size information could be extracted.
       **/
      template <class LF, class... Args,
        IsLocalFunction<LF,Element> = true,
        class R = Range<LF,LocalDomain> >
      Function (LF&& localFct, std::string name, Args&&... args)
        : Function(std::forward<LF>(localFct), std::move(name), sizeOf<R>(),
                   std::forward<Args>(args)...)
      {}

      /// (4) Construct from a Vtk::Function
      template <class... Args>
      explicit Function (Function<GridView> const& fct, Args&&... args)
        : Function(fct.localFct_,
                   getArg<std::string, char const*>(args..., fct.name_),
                   getArg<int,unsigned int,long,unsigned long,std::vector<int>>(args..., fct.components_),
                   getArg<Vtk::RangeTypes>(args..., fct.rangeType_),
                   getArg<Vtk::DataTypes>(args..., fct.dataType_))
      {}

      /// (5) Construct from a GridFunction
      /**
      * \param fct   A Grid(View)-function, providing a `localFunction(fct)`
      * \param name  The name to use as identification in the VTK file
      *
      * Forwards all other arguments to the constructor (1) or (2).
      *
      * NOTE: Stores the localFunction(fct) by value.
      */
      template <class GF, class... Args,
        disableCopyMove<Function, GF> = 0,
        IsGridFunction<GF> = true>
      Function (GF&& fct, std::string name, Args&&... args)
        : Function(localFunction(std::forward<GF>(fct)), std::move(name), std::forward<Args>(args)...)
      {}

      /// (6) Constructor that forwards the number of components and data type to the other constructor
      template <class F>
      Function (F&& fct, Vtk::FieldInfo info, ...)
        : Function(std::forward<F>(fct), info.name(), info.size(), info.rangeType(), info.dataType())
      {}

      /// (7) Automatically extract name and num components from GridFunction if available
      template <class F, class... Args,
        disableCopyMove<Function, F> = 0,
        IsGridFunction<F> = true,
        class = decltype(std::declval<F>().name()),
        class = decltype(std::declval<F>().numComponents()),
        class = decltype(std::declval<F>().dataType()) >
      explicit Function (F&& fct, ...)
        : Function(localFunction(std::forward<F>(fct)), fct.name(), fct.numComponents(),
                   Vtk::RangeTypes::UNSPECIFIED, fct.dataType())
      {}

      /// (8) Construct from legacy VTKFunction
      /**
      * \param fct  The Dune::VTKFunction to wrap
      **/
      explicit Function (std::shared_ptr<VTKFunction<GridView> const> const& fct, ...)
        : localFct_(fct)
        , name_(fct->name())
      {
        setComponents(fct->ncomps());
        setDataType(dataTypeOf(fct->precision()));
        setRangeType(rangeTypeOf(fct->ncomps()));
      }

      /// (9) Default constructor. After construction, the function is an an invalid state.
      Function () = default;

      /// Create a LocalFunction
      friend Vtk::LocalFunction<GridView> localFunction (Function const& self)
      {
        return self.localFct_;
      }

      /// Return a name associated with the function
      std::string const& name () const
      {
        return name_;
      }

      /// Set the function name
      void setName (std::string name)
      {
        name_ = std::move(name);
      }

      /// Return the number of components of the Range as it is written to the file
      int numComponents () const
      {
        return rangeType_ == Vtk::RangeTypes::SCALAR ? 1 :
               rangeType_ == Vtk::RangeTypes::VECTOR ? 3 :
               rangeType_ == Vtk::RangeTypes::TENSOR ? 9 : int(components_.size());
      }

      /// Set the components of the Range to visualize
      void setComponents (std::vector<int> components)
      {
        components_ = components;
        localFct_.setComponents(components_);
      }

      /// Set the number of components of the Range and generate component range [0...ncomps)
      void setComponents (int ncomps)
      {
        setComponents(allComponents(ncomps));
      }

      /// Return the VTK Datatype associated with the functions range type
      Vtk::DataTypes dataType () const
      {
        return dataType_;
      }

      /// Set the data-type for the components
      void setDataType (Vtk::DataTypes type)
      {
        dataType_ = type;
      }

      /// The category of the range, SCALAR, VECTOR, TENSOR, or UNSPECIFIED
      Vtk::RangeTypes rangeType () const
      {
        return rangeType_;
      }

      /// Set the category of the range, SCALAR, VECTOR, TENSOR, or UNSPECIFIED
      void setRangeType (Vtk::RangeTypes type, std::size_t ncomp = 1)
      {
        rangeType_ = type;
        if (type == Vtk::RangeTypes::AUTO)
          rangeType_ = rangeTypeOf(ncomp);
      }

      /// Set all the parameters from a FieldInfo object
      void setFieldInfo (Vtk::FieldInfo info)
      {
        setName(info.name());
        setComponents(info.size());
        setRangeType(info.rangeType());
        setDataType(info.dataType());
      }

    private:
      Vtk::LocalFunction<GridView> localFct_;
      std::string name_;
      std::vector<int> components_;
      Vtk::DataTypes dataType_ = Vtk::DataTypes::FLOAT64;
      Vtk::RangeTypes rangeType_ = Vtk::RangeTypes::UNSPECIFIED;
    };

  } // end namespace Vtk
} // end namespace Dune
