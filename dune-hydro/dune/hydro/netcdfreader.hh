#ifndef Udune_hydro_netcdfreader_HH
#define Udune_hydro_netcdfreader_HH

#include <exception>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

// netCDF include
#ifdef NetCDF_FOUND
#include <netcdf.h>
#endif

class NetCDFFileException : public std::exception
{
public:
  virtual const char* what() const throw()
  {
    return message.c_str();
  }

  NetCDFFileException (std::string m) : message(m) {}
  
private:
  std::string message;
};


class NetCDFFile
{
  std::string name; // name of the file (without the path)
  int verbose;
  int ncid; // the file id in the netcdf library
  int formatid; // which format is it?
  int ndims; // number of dimensions
  int nvars; // number of variables
  int ngatts; // number of global attributes
  int unlimdimid; // id of unlimited direction
  std::vector<std::string> name_of_dim;
  std::vector<size_t> length_of_dim;
  std::vector<std::string> name_of_var;
  std::vector<std::vector<int>> dimensions_of_var;
  std::vector<size_t> number_of_attributes;
  std::vector<nc_type> type_of_var;
  std::vector<std::vector<std::string>> var_attribute_name;
  std::vector<std::vector<nc_type>> var_attribute_type;
  std::vector<std::vector<size_t>> var_attribute_length;
  std::vector<std::string> global_attribute_name;
  std::vector<nc_type> global_attribute_type;
  std::vector<size_t> global_attribute_length;

public:

  NetCDFFile (std::string filename, std::string path="",int _verbose=2)
    : name(filename), verbose(_verbose)
  {
    /* Open the file. */
    int retval;
    std::string fullname;
    if (path.size()==0)
      fullname = filename;
    else
      fullname = path+filename;
    if ((retval = nc_open(fullname.c_str(), NC_NOWRITE, &ncid)))
      throw NetCDFFileException(std::string(nc_strerror(retval)));


    if (retval = nc_inq_format(ncid,&formatid))
      throw NetCDFFileException(std::string(nc_strerror(retval)));
    if (verbose>0)
      std::cout << "opened netcdf file " << filename << " with format " << format_to_string(formatid) << std::endl;

    /* There are a number of inquiry functions in netCDF which can be
       used to learn about an unknown netCDF file. NC_INQ tells how
       many netCDF variables, dimensions, and global attributes are in
       the file; also the dimension id of the unlimited dimension, if
       there is one. */
    if ((retval = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid)))
      throw NetCDFFileException(std::string(nc_strerror(retval)));


    // analyse the dimensions
    if (verbose>0)
      std::cout << ndims << " dimensions: " << std::endl;
    for (int dimid=0; dimid<ndims; dimid++)
      {
        size_t len;
        char name[NC_MAX_NAME];
        if (retval=nc_inq_dim(ncid,dimid,name,&len))
          throw NetCDFFileException(std::string(nc_strerror(retval)));
        name_of_dim.push_back(std::string(name));
        length_of_dim.push_back(len);
        if (verbose>1)
          std::cout << "  dim " << dimid << " is " << name << " with length " << len << std::endl;
      }
    if (verbose>1)
        std::cout << "  unlimited dimension is " << unlimdimid << std::endl;
    
    // analyse the variables
    if (verbose>0)
      std::cout << nvars << " variables: " << std::endl;
    if (nvars>0)
      {
        var_attribute_name.resize(nvars);
        var_attribute_type.resize(nvars);
        var_attribute_length.resize(nvars);
      }
    for (int varid=0; varid<nvars; varid++)
      {
        char name[NC_MAX_NAME];
        nc_type xtype;
        int ndims;
        int dimids[NC_MAX_VAR_DIMS];
        int natts;
        if (retval=nc_inq_var(ncid,varid,name,&xtype,&ndims,dimids,&natts))
          throw NetCDFFileException(std::string(nc_strerror(retval)));
        name_of_var.push_back(std::string(name));
        type_of_var.push_back(xtype);
        number_of_attributes.push_back(natts);
        std::vector<int> vardims;
        for (int i=0; i<ndims; i++) vardims.push_back(dimids[i]);
        dimensions_of_var.push_back(vardims);
        if (verbose>1)
          {
            std::cout << "  var " << varid << " of type " << nctype_to_string(xtype) << " is " << name << " with " << natts << " attributes and uses dimensions";
            for (int i=0; i<ndims; i++) std::cout << " " << dimids[i];
            std::cout << std::endl;
          }
        // analyse attributes of this variable
        for (int attid=0; attid<natts; attid++)
          {
            char name[NC_MAX_NAME]; name[0]=0;
            if (retval=nc_inq_attname(ncid,varid,attid,name))
              throw NetCDFFileException(std::string(nc_strerror(retval)));
            var_attribute_name[varid].push_back(std::string(name));
            nc_type xtype;
            size_t len;
            if (retval=nc_inq_att(ncid,varid,name,&xtype,&len))
              throw NetCDFFileException(std::string(nc_strerror(retval)));
            if (verbose>1)
              std::cout << "    attribute " << attid << " is " << name << " with type " << nctype_to_string(xtype) << " and length " << len << std::endl;
            var_attribute_type[varid].push_back(xtype);
            var_attribute_length[varid].push_back(len);
          }
      }

    // analyse global attributes
    if (verbose>0)
      std::cout << ngatts << " global attributes: " << std::endl;
    for (int attid=0; attid<ngatts; attid++)
      {
        char name[NC_MAX_NAME]; name[0]=0;
        if (retval=nc_inq_attname(ncid,NC_GLOBAL,attid,name))
          throw NetCDFFileException(std::string(nc_strerror(retval)));
        global_attribute_name.push_back(std::string(name));
        nc_type xtype;
        size_t len;
        if (retval=nc_inq_att(ncid,NC_GLOBAL,name,&xtype,&len))
          throw NetCDFFileException(std::string(nc_strerror(retval)));
        global_attribute_type.push_back(xtype);
        global_attribute_length.push_back(len);
        if (verbose>1)
          std::cout << "  global attribute " << attid << " is " << name << " with type " << nctype_to_string(xtype) << " and length " << len << std::endl;
      }
  }

  //! \brief get file format id
  int file_format () const
  {
    return formatid;
  }

  //! \brief  number of available dimensions in the data set
  int dimensions () const
  {
    return ndims;
  }

  //! \brief get name of a dimension
  std::string dimension_name (int dimid) const
  {
    return name_of_dim[dimid];
  }

  //! \brief get length of the dimension
  size_t dimension_length (int i) const
  {
    return length_of_dim[i];
  }

  //! \brief get id if the unlimited dimension; -1 if there is none
  int unlimited_dimension_id () const
  {
    return unlimdimid;
  }
  
  //! \brief number of variables in the data set
  int variables () const
  {
    return nvars;
  }

  //! \brief get name of the variable
  std::string variable_name (int i) const
  {
    return name_of_var[i];
  }

  //! \brief get type of the variable
  nc_type variable_type (int i) const
  {
    return type_of_var[i];
  }

  //! \brief get number of dimensions of the variable
  size_t variable_dimensions (int varid) const
  {
    return dimensions_of_var[varid].size();
  }

  //! \brief get dimension id of a dimension of a variable 
  int variable_dimension (int varid, int dimid) const
  {
    return dimensions_of_var[varid][dimid];
  }


  //! \brief number of global attributes
  int attributes () const
  {
    return ngatts;
  }

  //! \brief number of attributes of a variable
  int attributes (int varid) const
  {
    return var_attribute_name[varid].size();
  }

  //! \brief type of global attribute
  int attribute_type (int attid) const
  {
    return global_attribute_type[attid];
  }

  //! \brief type of variable attribute
  int attribute_type (int varid, int attid) const
  {
    return var_attribute_type[varid][attid];
  }

  //! \brief length of global attribute
  int attribute_length (int attid) const
  {
    return global_attribute_length[attid];
  }

  //! \brief type of variable attribute
  int attribute_length (int varid, int attid) const
  {
    return var_attribute_length[varid][attid];
  }

  //! \brief get global attribute of type string
  std::string get_string_attribute (int attid)  const
  {
    // check if this attribute is really a string
    if (attid>=global_attribute_type.size())
      throw NetCDFFileException("get_global_attribute: invalid attid");
    if (global_attribute_type[attid]!=NC_CHAR)
      throw NetCDFFileException("get_global_attribute: data type is not char");

    char *p = new char[global_attribute_length[attid]+5];
    auto retval = nc_get_att_text(ncid,NC_GLOBAL,global_attribute_name[attid].c_str(),p);
    p[global_attribute_length[attid]] = 0;
    std::string s(p);
    delete[] p;
    return s;
  }

  // get variable attribute of type string
  std::string get_string_attribute (int varid, int attid) const
  {
    // check if this attribute is really a string
    if (varid>=nvars)
      throw NetCDFFileException("get_global_attribute: invalid varid");
    if (attid>=number_of_attributes[varid])
      throw NetCDFFileException("get_global_attribute: invalid attid");
    if (var_attribute_type[varid][attid]!=NC_CHAR)
      throw NetCDFFileException("get_global_attribute: data type is not char");

    char *p = new char[var_attribute_length[varid][attid]+5];
    auto retval = nc_get_att_text(ncid,varid,var_attribute_name[varid][attid].c_str(),p);
    p[var_attribute_length[varid][attid]] = 0;
    std::string s(p);
    delete[] p;
    return s;
  }

  // read variable attribute of size 1
  template<typename T>
  T get_attribute (int varid, int attid) const
  {
    if (varid>=nvars) // check varid
      throw NetCDFFileException("invalid varid");
    if (attid>=number_of_attributes[varid]) // check attid
      throw NetCDFFileException("invalid attid");
    if (var_attribute_length[varid][attid]!=1) // check size is 1
      throw NetCDFFileException("works only for size 1");
    if (var_attribute_type[varid][attid]==NC_CHAR)
      throw NetCDFFileException("use get_string_attribute for char");
    if (var_attribute_type[varid][attid]==NC_BYTE && std::is_same<T,char>::value==false)
      throw NetCDFFileException("data type char was expected");
    if (var_attribute_type[varid][attid]==NC_SHORT && std::is_same<T,short>::value==false)
      throw NetCDFFileException("data type short was expected");
    if (var_attribute_type[varid][attid]==NC_INT && std::is_same<T,int>::value==false)
      throw NetCDFFileException("data type int was expected");
    if (var_attribute_type[varid][attid]==NC_LONG && std::is_same<T,long>::value==false)
      throw NetCDFFileException("data type long was expected");
    if (var_attribute_type[varid][attid]==NC_INT64 && std::is_same<T,long>::value==false)
      throw NetCDFFileException("data type long was expected");
    if (var_attribute_type[varid][attid]==NC_FLOAT && std::is_same<T,float>::value==false)
      throw NetCDFFileException("data type float was expected");
    if (var_attribute_type[varid][attid]==NC_DOUBLE && std::is_same<T,double>::value==false)
      throw NetCDFFileException("data type double was expected");
    if (var_attribute_type[varid][attid]==NC_UBYTE && std::is_same<T,unsigned char>::value==false)
      throw NetCDFFileException("data type unsigned char was expected");
    if (var_attribute_type[varid][attid]==NC_USHORT && std::is_same<T,unsigned short>::value==false)
      throw NetCDFFileException("data type unsigned short was expected");
    if (var_attribute_type[varid][attid]==NC_UINT && std::is_same<T,unsigned int>::value==false)
      throw NetCDFFileException("data type unsigned int was expected");
    if (var_attribute_type[varid][attid]==NC_UINT64 && std::is_same<T,unsigned long>::value==false)
      throw NetCDFFileException("data type unsigned long was expected");

    T data;
    nc_get_att(ncid,varid,var_attribute_name[varid][attid].c_str(),(void *) &data);
    return data;
  }

  // read global attribute of size 1
  template<typename T>
  T get_attribute (int attid) const
  {
    if (attid>=global_attribute_type.size()) // check attid
      throw NetCDFFileException("invalid attid");
    if (global_attribute_length[attid]!=1) // check size is 1
      throw NetCDFFileException("works only for size 1");
    if (global_attribute_type[attid]==NC_CHAR)
      throw NetCDFFileException("use get_string_attribute for char");
    if (global_attribute_type[attid]==NC_BYTE && std::is_same<T,char>::value==false)
      throw NetCDFFileException("data type char was expected");
    if (global_attribute_type[attid]==NC_SHORT && std::is_same<T,short>::value==false)
      throw NetCDFFileException("data type short was expected");
    if (global_attribute_type[attid]==NC_INT && std::is_same<T,int>::value==false)
      throw NetCDFFileException("data type int was expected");
    if (global_attribute_type[attid]==NC_LONG && std::is_same<T,long>::value==false)
      throw NetCDFFileException("data type long was expected");
    if (global_attribute_type[attid]==NC_INT64 && std::is_same<T,long>::value==false)
      throw NetCDFFileException("data type long was expected");
    if (global_attribute_type[attid]==NC_FLOAT && std::is_same<T,float>::value==false)
      throw NetCDFFileException("data type float was expected");
    if (global_attribute_type[attid]==NC_DOUBLE && std::is_same<T,double>::value==false)
      throw NetCDFFileException("data type double was expected");
    if (global_attribute_type[attid]==NC_UBYTE && std::is_same<T,unsigned char>::value==false)
      throw NetCDFFileException("data type unsigned char was expected");
    if (global_attribute_type[attid]==NC_USHORT && std::is_same<T,unsigned short>::value==false)
      throw NetCDFFileException("data type unsigned short was expected");
    if (global_attribute_type[attid]==NC_UINT && std::is_same<T,unsigned int>::value==false)
      throw NetCDFFileException("data type unsigned int was expected");
    if (global_attribute_type[attid]==NC_UINT64 && std::is_same<T,unsigned long>::value==false)
      throw NetCDFFileException("data type unsigned long was expected");

    T data;
    nc_get_att(ncid,NC_GLOBAL,global_attribute_name[attid].c_str(),(void *) &data);
    return data;
  }

  // read variable attribute vector
  template<typename T>
  std::vector<T> get_attribute_vector (int varid, int attid) const
  {
    if (varid>=nvars) // check varid
      throw NetCDFFileException("invalid varid");
    if (attid>=number_of_attributes[varid]) // check attid
      throw NetCDFFileException("invalid attid");
    if (var_attribute_type[varid][attid]==NC_CHAR)
      throw NetCDFFileException("use get_string_attribute for char");
    if (var_attribute_type[varid][attid]==NC_BYTE && std::is_same<T,char>::value==false)
      throw NetCDFFileException("data type char was expected");
    if (var_attribute_type[varid][attid]==NC_SHORT && std::is_same<T,short>::value==false)
      throw NetCDFFileException("data type short was expected");
    if (var_attribute_type[varid][attid]==NC_INT && std::is_same<T,int>::value==false)
      throw NetCDFFileException("data type int was expected");
    if (var_attribute_type[varid][attid]==NC_LONG && std::is_same<T,long>::value==false)
      throw NetCDFFileException("data type long was expected");
    if (var_attribute_type[varid][attid]==NC_INT64 && std::is_same<T,long>::value==false)
      throw NetCDFFileException("data type long was expected");
    if (var_attribute_type[varid][attid]==NC_FLOAT && std::is_same<T,float>::value==false)
      throw NetCDFFileException("data type float was expected");
    if (var_attribute_type[varid][attid]==NC_DOUBLE && std::is_same<T,double>::value==false)
      throw NetCDFFileException("data type double was expected");
    if (var_attribute_type[varid][attid]==NC_UBYTE && std::is_same<T,unsigned char>::value==false)
      throw NetCDFFileException("data type unsigned char was expected");
    if (var_attribute_type[varid][attid]==NC_USHORT && std::is_same<T,unsigned short>::value==false)
      throw NetCDFFileException("data type unsigned short was expected");
    if (var_attribute_type[varid][attid]==NC_UINT && std::is_same<T,unsigned int>::value==false)
      throw NetCDFFileException("data type unsigned int was expected");
    if (var_attribute_type[varid][attid]==NC_UINT64 && std::is_same<T,unsigned long>::value==false)
      throw NetCDFFileException("data type unsigned long was expected");

    std::vector<T> v(var_attribute_length[varid][attid]);
    nc_get_att(ncid,varid,var_attribute_name[varid][attid].c_str(),(void *) v.data());
    return v;
  }

  // read global attribute of size 1
  template<typename T>
  T get_attribute_vector (int attid) const
  {
    if (attid>=global_attribute_type.size()) // check attid
      throw NetCDFFileException("invalid attid");
    if (global_attribute_type[attid]==NC_CHAR)
      throw NetCDFFileException("use get_string_attribute for char");
    if (global_attribute_type[attid]==NC_BYTE && std::is_same<T,char>::value==false)
      throw NetCDFFileException("data type char was expected");
    if (global_attribute_type[attid]==NC_SHORT && std::is_same<T,short>::value==false)
      throw NetCDFFileException("data type short was expected");
    if (global_attribute_type[attid]==NC_INT && std::is_same<T,int>::value==false)
      throw NetCDFFileException("data type int was expected");
    if (global_attribute_type[attid]==NC_LONG && std::is_same<T,long>::value==false)
      throw NetCDFFileException("data type long was expected");
    if (global_attribute_type[attid]==NC_INT64 && std::is_same<T,long>::value==false)
      throw NetCDFFileException("data type long was expected");
    if (global_attribute_type[attid]==NC_FLOAT && std::is_same<T,float>::value==false)
      throw NetCDFFileException("data type float was expected");
    if (global_attribute_type[attid]==NC_DOUBLE && std::is_same<T,double>::value==false)
      throw NetCDFFileException("data type double was expected");
    if (global_attribute_type[attid]==NC_UBYTE && std::is_same<T,unsigned char>::value==false)
      throw NetCDFFileException("data type unsigned char was expected");
    if (global_attribute_type[attid]==NC_USHORT && std::is_same<T,unsigned short>::value==false)
      throw NetCDFFileException("data type unsigned short was expected");
    if (global_attribute_type[attid]==NC_UINT && std::is_same<T,unsigned int>::value==false)
      throw NetCDFFileException("data type unsigned int was expected");
    if (global_attribute_type[attid]==NC_UINT64 && std::is_same<T,unsigned long>::value==false)
      throw NetCDFFileException("data type unsigned long was expected");

    std::vector<T> v(global_attribute_length[attid]);
    nc_get_att(ncid,NC_GLOBAL,global_attribute_name[attid].c_str(),(void *) v.data());
    return v;
  }

  //! \brief read complete data of one variable
  template<typename T>
  std::vector<T> get_variable (int varid) const
  {
    if (varid>=nvars) // check attid
      throw NetCDFFileException("invalid varid");
    if (type_of_var[varid]==NC_CHAR)
      throw NetCDFFileException("use get_string_attribute for char");
    if (type_of_var[varid]==NC_BYTE && std::is_same<T,char>::value==false)
      throw NetCDFFileException("data type char was expected");
    if (type_of_var[varid]==NC_SHORT && std::is_same<T,short>::value==false)
      throw NetCDFFileException("data type short was expected");
    if (type_of_var[varid]==NC_INT && std::is_same<T,int>::value==false)
      throw NetCDFFileException("data type int was expected");
    if (type_of_var[varid]==NC_LONG && std::is_same<T,long>::value==false)
      throw NetCDFFileException("data type long was expected");
    if (type_of_var[varid]==NC_INT64 && std::is_same<T,long>::value==false)
      throw NetCDFFileException("data type long was expected");
    if (type_of_var[varid]==NC_FLOAT && std::is_same<T,float>::value==false)
      throw NetCDFFileException("data type float was expected");
    if (type_of_var[varid]==NC_DOUBLE && std::is_same<T,double>::value==false)
      throw NetCDFFileException("data type double was expected");
    if (type_of_var[varid]==NC_UBYTE && std::is_same<T,unsigned char>::value==false)
      throw NetCDFFileException("data type unsigned char was expected");
    if (type_of_var[varid]==NC_USHORT && std::is_same<T,unsigned short>::value==false)
      throw NetCDFFileException("data type unsigned short was expected");
    if (type_of_var[varid]==NC_UINT && std::is_same<T,unsigned int>::value==false)
      throw NetCDFFileException("data type unsigned int was expected");
    if (type_of_var[varid]==NC_UINT64 && std::is_same<T,unsigned long>::value==false)
      throw NetCDFFileException("data type unsigned long was expected");

    size_t n=1;
    for (int i=0; i<dimensions_of_var[varid].size(); i++)
      n *= length_of_dim[dimensions_of_var[varid][i]];
    if (verbose>1) std::cout << "total size is " << n << std::endl;
    std::vector<T> v(n);
    nc_get_var(ncid,varid,(void *) v.data());
    return v;
  }

  //! \brief read data w.r.t. to one value of the unlimited id of a variable
  template<typename T>
  std::vector<T> get_variable (int varid, int unlim_i) const
  {
    if (varid>=nvars) // check attid
      throw NetCDFFileException("invalid varid");
    if (unlimdimid<0) // check if there is an unlimited dimension
      throw NetCDFFileException("no unlimited dimension");
    if (dimensions_of_var[varid][0]!=unlimdimid)
      throw NetCDFFileException("unlimited dimension must be the slowest dimension of the variable");
    if (type_of_var[varid]==NC_CHAR)
      throw NetCDFFileException("use get_string_attribute for char");
    if (type_of_var[varid]==NC_BYTE && std::is_same<T,char>::value==false)
      throw NetCDFFileException("data type char was expected");
    if (type_of_var[varid]==NC_SHORT && std::is_same<T,short>::value==false)
      throw NetCDFFileException("data type short was expected");
    if (type_of_var[varid]==NC_INT && std::is_same<T,int>::value==false)
      throw NetCDFFileException("data type int was expected");
    if (type_of_var[varid]==NC_LONG && std::is_same<T,long>::value==false)
      throw NetCDFFileException("data type long was expected");
    if (type_of_var[varid]==NC_INT64 && std::is_same<T,long>::value==false)
      throw NetCDFFileException("data type long was expected");
    if (type_of_var[varid]==NC_FLOAT && std::is_same<T,float>::value==false)
      throw NetCDFFileException("data type float was expected");
    if (type_of_var[varid]==NC_DOUBLE && std::is_same<T,double>::value==false)
      throw NetCDFFileException("data type double was expected");
    if (type_of_var[varid]==NC_UBYTE && std::is_same<T,unsigned char>::value==false)
      throw NetCDFFileException("data type unsigned char was expected");
    if (type_of_var[varid]==NC_USHORT && std::is_same<T,unsigned short>::value==false)
      throw NetCDFFileException("data type unsigned short was expected");
    if (type_of_var[varid]==NC_UINT && std::is_same<T,unsigned int>::value==false)
      throw NetCDFFileException("data type unsigned int was expected");
    if (type_of_var[varid]==NC_UINT64 && std::is_same<T,unsigned long>::value==false)
      throw NetCDFFileException("data type unsigned long was expected");

    size_t n=1; // the size of the buffer to be returned
    for (int i=1; i<dimensions_of_var[varid].size(); i++)
      n *= length_of_dim[dimensions_of_var[varid][i]];
    if (verbose>1) std::cout << "total size is " << n << std::endl;
    std::vector<T> v(n);
    std::vector<size_t> start(dimensions_of_var[varid].size(),0);
    start[0] = unlim_i;
    std::vector<size_t> count(dimensions_of_var[varid].size(),1);
    for (int i=1; i<dimensions_of_var[varid].size(); i++)
      count[i] = length_of_dim[dimensions_of_var[varid][i]];
    nc_get_vara(ncid,varid,start.data(),count.data(),(void *) v.data());
    return v;
  }
  
  ~NetCDFFile ()
  {
    /* Close the file. */
    nc_close(ncid);
  }

  std::string format_to_string (int f) const
  {
    if (f==NC_FORMAT_CLASSIC) return std::string("NC_FORMAT_CLASSIC");
    if (f==NC_FORMAT_64BIT_OFFSET) return std::string("NC_FORMAT_64BIT_OFFSET");
    if (f==NC_FORMAT_CDF5) return std::string("NC_FORMAT_CDF5");
    if (f==NC_FORMAT_NETCDF4) return std::string("NC_FORMAT_NETCDF4");
    if (f==NC_FORMAT_NETCDF4_CLASSIC) return std::string("NC_FORMAT_NETCDF4_CLASSIC");
    return std::string("NONE!");
  }

  std::string nctype_to_string (nc_type t) const
  {
    if (t==NC_NAT) return std::string("NC_NAT");
    if (t==NC_BYTE) return std::string("NC_BYTE");
    if (t==NC_CHAR) return std::string("NC_CHAR");
    if (t==NC_SHORT) return std::string("NC_SHORT");
    if (t==NC_INT) return std::string("NC_INT");
    if (t==NC_LONG ) return std::string("NC_LONG");
    if (t==NC_FLOAT ) return std::string("NC_FLOAT");
    if (t==NC_DOUBLE) return std::string("NC_DOUBLE");
    if (t==NC_UBYTE) return std::string("NC_UBYTE");
    if (t==NC_USHORT) return std::string("NC_USHORT");
    if (t==NC_UINT) return std::string("NC_UINT");
    if (t==NC_INT64) return std::string("NC_INT64");
    if (t==NC_UINT64) return std::string("NC_UINT64");
    if (t==NC_STRING) return std::string("NC_STRING");
    if (t==NC_MAX_ATOMIC_TYPE) return std::string("NC_MAX_ATOMIC_TYPE");
    return std::string("NONE!");
  }
};

template<typename T>
void print_vector_overview (const std::vector<T> c, int n=4)
{
  for (int i=0; i<std::min((size_t)n,c.size()); i++)
    std::cout << c[i] << " ";
  std::cout << "...";
  for (int i=std::max((size_t)0,c.size()-n); i<c.size(); i++)
    std::cout << " " << c[i];
  std::cout << std::endl;
}

#endif 
