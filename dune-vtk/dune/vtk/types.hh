#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include <dune/common/ftraits.hh>
#include <dune/common/typelist.hh>
#include <dune/common/version.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/vtk/utility/arguments.hh>
#include <dune/vtk/utility/errors.hh>

namespace Dune
{
  namespace Vtk
  {
    /// Type used for representing the output format
    enum class FormatTypes {
      ASCII      = 1<<0,
      BINARY     = 1<<1,
      COMPRESSED = 1<<2,
      APPENDED = BINARY | COMPRESSED
    };
    std::string to_string (Vtk::FormatTypes);
    inline auto formatTypesList = {FormatTypes::ASCII, FormatTypes::BINARY, FormatTypes::COMPRESSED, FormatTypes::APPENDED};

    /// Map the dune-grid OutputType to FormatTypes
    Vtk::FormatTypes formatTypeOf(Dune::VTK::OutputType);


    /// Type used to determine whether to limit output components to e.g. 3 (vector), or 9 (tensor)
    enum class RangeTypes {
      UNSPECIFIED,  //< The output components are not restricted
      AUTO,         //< Detect the category automatically from number of components
      SCALAR,       //< Use exactly 1 component
      VECTOR,       //< Use exactly 3 components
      TENSOR        //< Use exactly 9 components
    };
    std::string to_string (Vtk::RangeTypes);
    inline auto rangeTypesList = {RangeTypes::UNSPECIFIED, RangeTypes::AUTO, RangeTypes::SCALAR, RangeTypes::VECTOR, RangeTypes::TENSOR};

    // Map a dune-grid FieldInfo::Type to ValueTypes
    Vtk::RangeTypes rangeTypeOf (Dune::VTK::FieldInfo::Type);

    // Map a number of components to a corresponding value type
    Vtk::RangeTypes rangeTypeOf (int ncomps);


    enum class DataTypes {
      UNKNOWN = 0,
      INT8, UINT8,
      INT16, UINT16,
      INT32, UINT32,
      INT64, UINT64,
      FLOAT32 = 32,
      FLOAT64 = 64
    };
    std::string to_string (Vtk::DataTypes);
    inline auto dataTypesLists = {
      DataTypes::UNKNOWN,
      DataTypes::INT8,    DataTypes::UINT8,
      DataTypes::INT16,   DataTypes::UINT16,
      DataTypes::INT32,   DataTypes::UINT32,
      DataTypes::INT64,   DataTypes::UINT64,
      DataTypes::FLOAT32, DataTypes::FLOAT64
    };

    // Map a dune-grid Precision type to DataTypes
    Vtk::DataTypes dataTypeOf (Dune::VTK::Precision);

    // Map a string to DataTypes
    Vtk::DataTypes dataTypeOf (std::string);

    // Map the field_type of T to DataTypes
    template <class T>
    Vtk::DataTypes dataTypeOf ()
    {
      using F = typename FieldTraits<T>::field_type;
      if constexpr (std::is_same_v<F, std::int8_t>)   { return Vtk::DataTypes::INT8; }
      if constexpr (std::is_same_v<F, std::uint8_t>)  { return Vtk::DataTypes::UINT8; }
      if constexpr (std::is_same_v<F, std::int16_t>)  { return Vtk::DataTypes::INT16; }
      if constexpr (std::is_same_v<F, std::uint16_t>) { return Vtk::DataTypes::UINT16; }
      if constexpr (std::is_same_v<F, std::int32_t>)  { return Vtk::DataTypes::INT32; }
      if constexpr (std::is_same_v<F, std::uint32_t>) { return Vtk::DataTypes::UINT32; }
      if constexpr (std::is_same_v<F, std::int64_t>)  { return Vtk::DataTypes::INT64; }
      if constexpr (std::is_same_v<F, std::uint64_t>) { return Vtk::DataTypes::UINT64; }
      if constexpr (std::is_same_v<F, float>)         { return Vtk::DataTypes::FLOAT32; }
      if constexpr (std::is_same_v<F, double>)        { return Vtk::DataTypes::FLOAT64; }
      if constexpr (std::is_same_v<F, long double>)   { return Vtk::DataTypes::FLOAT64; }
      return Vtk::DataTypes::UNKNOWN;
    }

    template <class> struct NoConstraint : std::true_type {};

    /// Map a given enum DataType to a type passed to Caller as \ref MetaType
    template <template <class> class C = NoConstraint, class Caller>
    void mapDataTypes (Vtk::DataTypes t, Caller caller)
    {
      switch (t) {
        case DataTypes::INT8:    if constexpr(C<std::int8_t>::value)   caller(MetaType<std::int8_t>{});   break;
        case DataTypes::UINT8:   if constexpr(C<std::uint8_t>::value)  caller(MetaType<std::uint8_t>{});  break;
        case DataTypes::INT16:   if constexpr(C<std::int16_t>::value)  caller(MetaType<std::int16_t>{});  break;
        case DataTypes::UINT16:  if constexpr(C<std::uint16_t>::value) caller(MetaType<std::uint16_t>{}); break;
        case DataTypes::INT32:   if constexpr(C<std::int32_t>::value)  caller(MetaType<std::int32_t>{});  break;
        case DataTypes::UINT32:  if constexpr(C<std::uint32_t>::value) caller(MetaType<std::uint32_t>{}); break;
        case DataTypes::INT64:   if constexpr(C<std::int64_t>::value)  caller(MetaType<std::int64_t>{});  break;
        case DataTypes::UINT64:  if constexpr(C<std::uint64_t>::value) caller(MetaType<std::uint64_t>{}); break;
        case DataTypes::FLOAT32: if constexpr(C<float>::value)         caller(MetaType<float>{});         break;
        case DataTypes::FLOAT64: if constexpr(C<double>::value)        caller(MetaType<double>{});        break;
        default:
          VTK_ASSERT_MSG(false, "Unsupported type " + to_string(t));
          break;
      }
    }

    /// Map two DataTypes as type parameters to the Caller
    template <template <class> class Constraint1 = NoConstraint,
              template <class> class Constraint2 = NoConstraint,
              class Caller>
    void mapDataTypes (Vtk::DataTypes t1, Vtk::DataTypes t2, Caller caller)
    {
      mapDataTypes<Constraint1>(t1, [&](auto type1) {
        mapDataTypes<Constraint2>(t2, [&](auto type2) {
          caller(type1, type2);
        });
      });
    }

    /// Map three DataTypes as type parameters to the Caller
    template <template <class> class Constraint1 = NoConstraint,
              template <class> class Constraint2 = NoConstraint,
              template <class> class Constraint3 = NoConstraint,
              class Caller>
    void mapDataTypes (Vtk::DataTypes t1, Vtk::DataTypes t2, Vtk::DataTypes t3, Caller caller)
    {
      mapDataTypes<Constraint1>(t1, [&](auto type1) {
        mapDataTypes<Constraint2>(t2, [&](auto type2) {
          mapDataTypes<Constraint3>(t3, [&](auto type3) {
            caller(type1, type2, type3);
          });
        });
      });
    }


    enum class CompressorTypes {
      NONE = 0,
      ZLIB,
      LZ4,
      LZMA
    };
    std::string to_string (CompressorTypes);


    /// Mapping of Dune geometry types to VTK cell types
    struct CellType
    {
      enum Parametrization {
        LINEAR      = 1,
        QUADRATIC   = 2,
        LAGRANGE    = 3
      };

      enum Type : std::uint8_t {
        // Linear VTK cell types
        VERTEX                = 1,
        /* POLY_VERTEX        = 2, // not supported */
        LINE                  = 3,
        /* POLY_LINE          = 4, // not supported */
        TRIANGLE              = 5,
        /* TRIANGLE_STRIP     = 6, // not supported */
        POLYGON               = 7,
        /* PIXEL              = 8, // not supported */
        QUAD                  = 9,
        TETRA                 = 10,
        /* VOXEL              = 11, // not supported */
        HEXAHEDRON            = 12,
        WEDGE                 = 13,
        PYRAMID               = 14,
        // Quadratic VTK cell types
        QUADRATIC_EDGE        = 21,
        QUADRATIC_TRIANGLE    = 22,
        QUADRATIC_QUAD        = 23,
        QUADRATIC_TETRA       = 24,
        QUADRATIC_HEXAHEDRON  = 25,
        // Arbitrary order Lagrange elements
        LAGRANGE_CURVE        = 68,
        LAGRANGE_TRIANGLE     = 69,
        LAGRANGE_QUADRILATERAL= 70,
        LAGRANGE_TETRAHEDRON  = 71,
        LAGRANGE_HEXAHEDRON   = 72,
        LAGRANGE_WEDGE        = 73,
        LAGRANGE_PYRAMID      = 74,
      };

    public:
      CellType (GeometryType const& t, Parametrization = LINEAR);

      /// Return VTK Cell type
      std::uint8_t type () const
      {
        return type_;
      }

      /// Return a permutation of Dune elemenr vertices to conform to VTK element numbering
      int permutation (int idx) const
      {
        return permutation_[idx];
      }

      bool noPermutation () const
      {
        return noPermutation_;
      }

    private:
      std::uint8_t type_;
      std::vector<int> permutation_;
      bool noPermutation_ = true;
    };
    GeometryType to_geometry (std::uint8_t);


    class FieldInfo
    {
    public:
      template <class... Args>
      explicit FieldInfo (std::string name, Args... args)
        : name_(std::move(name))
        , ncomps_(getArg<int,unsigned int,long,unsigned long>(args..., 1))
        , rangeType_(getArg<Vtk::RangeTypes>(args..., Vtk::RangeTypes::AUTO))
        , dataType_(getArg<Vtk::DataTypes>(args..., Vtk::DataTypes::FLOAT32))
      {
        if (rangeType_ == Vtk::RangeTypes::AUTO)
          rangeType_ = rangeTypeOf(ncomps_);
      }

      // Construct from dune-grid FieldInfo
      FieldInfo (Dune::VTK::FieldInfo info)
        : FieldInfo(info.name(), info.size(), rangeTypeOf(info.type()), dataTypeOf(info.precision()))
      {}

      /// The name of the data field
      std::string const& name () const
      {
        return name_;
      }

      /// The number of components in the data field.
      int size () const
      {
        return ncomps_;
      }

      /// Return the category of the stored range
      Vtk::RangeTypes rangeType () const
      {
        return rangeType_;
      }

      /// Return the data tpe of the data field.
      Vtk::DataTypes dataType () const
      {
        return dataType_;
      }

    private:
      std::string name_;
      int ncomps_;
      Vtk::RangeTypes rangeType_;
      Vtk::DataTypes dataType_;
    };

  } // end namespace Vtk
} // end namespace Dune
