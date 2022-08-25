#ifndef Udune_hydro_geotiffreader_HH
#define Udune_hydro_geotiffreader_HH

#include <exception>
#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>

//GDAL includes. We require that GDAL is found
#ifdef GDAL_FOUND
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#endif

// add your classes here
//********************************************************************************
// GDAL stuff. to be put into a seperate file later ...
//********************************************************************************

/** list all tiles needed for a certain region
    The arguments are given in latitude and longitude:
    long_lower  longitude of most west point, west of Greenwich is negative, east is positive
    lat_lower   latitude of most south point, south of equator is negative, north of equator is positive
    long_upper  longitude of most east point, west of Greenwich is negative, east is positive
    lat_upper   latitude of most north point, south of equator is negative, north of equator is positive
 */
void listTiles (double long_lower, double lat_lower, double long_upper, double lat_upper)
{
  // tiles have the format [n|s]YYY_[e,w]XXX_1arc_v3.tif
  int xint_lower = floor(long_lower);
  int yint_lower = floor(lat_lower);
  int xint_upper = ceil(long_upper);
  int yint_upper = ceil(lat_upper);
  for (int y=yint_lower; y<yint_upper; y++)
    for (int x=xint_lower; x<xint_upper; x++)
      {
        std::stringstream latstring;
        if (y>=0) latstring << "n" << y; else latstring << "s" << std::abs(y);      
        std::stringstream longstring;
        if (x>=0) longstring << "e" << x; else longstring << "w" << std::abs(x);
        std::stringstream tilename;
        tilename << latstring.str() << "_" << longstring.str() << "_1arc_v3.tif"; 
        std::cout << "need tile " << x << "," << y << " to " << x+1 << "," << y+1 << "  " << tilename.str() << std::endl;
      }
}

class GeoTIFFReaderException : public std::exception
{
public:
  virtual const char* what() const throw()
  {
    return message.c_str();
  }

  GeoTIFFReaderException (std::string m) : message(m) {}
  
private:
  std::string message;
};

/** A class that opens a GeoTIFFImage and provides all the info without reading data 
 */
class GeoTIFFInfo {
public:

  GeoTIFFInfo (std::string filename, std::string path="", int verbose=0)
    : image_name(filename)
  {
    GDALDataset  *poDataset;
    std::string fullname;
    if (path.size()==0)
      fullname = filename;
    else
      fullname = path+filename;
    poDataset = (GDALDataset *) GDALOpen( fullname.c_str(), GA_ReadOnly );
    if ( poDataset == NULL )
      throw GeoTIFFReaderException("GDALOpen() failed");

    // read size
    xsize = poDataset->GetRasterXSize();
    ysize = poDataset->GetRasterYSize();
    count = poDataset->GetRasterCount();

    // read transformation to lat long coordinates
    double GT[6];
    if ( poDataset->GetGeoTransform( GT ) != CE_None )
      throw GeoTIFFReaderException("GetGeoTransform() failed");
    oX = GT[0];  oY = GT[3];
    AXX = GT[1]; AXY = GT[2];
    AYX = GT[4]; AYY = GT[5];
    cX = oX + 0.5*AXX; // half pixel offset since we report the center
    cY = oY + 0.5*AYY;
    firstX = std::min(cX,cX+(xsize-1)*AXX);
    firstY = std::min(cY,cY+(ysize-1)*AYY);
      
    if (verbose>0)
      {
	double extendx = xsize*AXX;
	double extendy = ysize*AYY;
	double originx = std::min(oX,oX+extendx);
	double originy = std::min(oY,oY+extendy);
	std::cout << "GeoTIFF file " << filename << " opened" << std::endl;
	std::cout << " size=" << xsize << "x" << ysize << " with " << count << " band(s)" << std::endl;
	std::cout << " origin=" << originx << "," << originy << " extend=" << std::abs(extendx) << "," << std::abs(extendy)
		  << " first pixel at " << firstX << "," << firstY << std::endl;
	if (std::abs(cX-round(cX))<1e-6 && std::abs(cY-round(cY))<1e-6)
	  std::cout << " vertex centered raster image" << std::endl;
	else
	  std::cout << " element centered raster image" << std::endl;	  
      }

    // read raster band
    GDALRasterBand* band = poDataset->GetRasterBand(1);
    dt = band->GetRasterDataType();
    _noDataValue = band->GetNoDataValue(&_hasNoDataValue);
    if (verbose>0)
      {
	std::cout << " DataType=" << dt << std::endl;
	if (_hasNoDataValue) std::cout << " no data value=" << _noDataValue << std::endl;
      }

    // read block information
    int nXBlockSize, nYBlockSize;
    band->GetBlockSize( &nXBlockSize, &nYBlockSize );
    if (verbose>0)
      std::cout << " nXBlockSize=" << nXBlockSize << " nYBlockSize=" << nYBlockSize << std::endl;
    int nXBlocks = (band->GetXSize() + nXBlockSize - 1) / nXBlockSize;
    int nYBlocks = (band->GetYSize() + nYBlockSize - 1) / nYBlockSize;
    if (verbose>0)
      std::cout << " nXBlocks=" << nXBlocks << " nYBlocks=" << nYBlocks << std::endl;
  }

  GDALDataType getGDALDataType ()
  {
    return dt;
  }
  
  //! apply affine GeoTransform to arbitrary vector data type
  template<typename U>
  U transform (const U& x) const
  {
    U y(x); // make a copy
    y[0] =  oX + AXX*x[0] + AXY*x[1];
    y[1] =  oY + AYX*x[0] + AYY*x[1];
    return y;
  }
  
  //! number of pixels in longitude direction
  int sizeLong () const
  {
    return xsize;
  }

  //! number of pixels in latitude direction
  int sizeLat () const
  {
    return ysize;
  }

  //! longitude of pixel *center*
  double longitude (int x, int y) const
  {
    return cX + AXX*x + AXY*y;
  }

  //! latitude of pixel *center*
  double latitude (int x, int y) const
  {
    return cY + AYX*x + AYY*y;
  }

  //! Pixel size in longitude direction
  double originLong () const
  {
    return firstX;
  }

  //! Pixel size in latitude direction
  double originLat () const
  {
    return firstY;
  }

  //! Pixel size in longitude direction
  double dLong () const
  {
    return AXX;
  }

  //! Pixel size in latitude direction
  double dLat () const
  {
    return AYY;
  }

  bool hasNoDataValue () const
  {
    return (_hasNoDataValue!=0);
  }

  double noDataValue () const
  {
    return _noDataValue;
  }
  
private:
  int xsize, ysize; // number of pixels in x and y
  int count; // number of raster bands in the data set
  double oX,oY;     // origin of image
  double cX,cY;     // origin of image offset by a half pixel
  double AXX, AXY, AYX, AYY; // transformation matrix
  std::string image_name; // name of the image
  double firstX, firstY; // position of first pixel in lower left corner (instead of up and left)!
  GDALDataType dt;
  int _hasNoDataValue;
  double _noDataValue;
};



/** A class that reads a GeoTIFFImage
    
    - holds the data in a buffer and lets you access the pixels
    - provides size of the image
    - provides positions of pixel centers in lat-long
    - Note: images are y-downwards!
 */
template<typename T=GInt16>
class GeoTIFFImage {
public:
  using value_type = T;

  GeoTIFFImage (std::string filename, std::string path="", int useband=1, int verbose=0)
    : image_name(filename)
  {
    GDALDataset  *poDataset;
    std::string fullname;
    if (path.size()==0)
      fullname = filename;
    else
      fullname = path+filename;
    poDataset = (GDALDataset *) GDALOpen( fullname.c_str(), GA_ReadOnly );
    if ( poDataset == NULL )
      throw GeoTIFFReaderException("GDALOpen() failed");

    // read size
    xsize = poDataset->GetRasterXSize();
    ysize = poDataset->GetRasterYSize();
    count = poDataset->GetRasterCount();

    // read transformation to lat long coordinates
    double GT[6];
    if ( poDataset->GetGeoTransform( GT ) != CE_None )
      throw GeoTIFFReaderException("GetGeoTransform() failed");
    oX = GT[0];  oY = GT[3];
    AXX = GT[1]; AXY = GT[2];
    AYX = GT[4]; AYY = GT[5];
    cX = oX + 0.5*AXX; // half pixel offset since we report the center
    cY = oY + 0.5*AYY;
    firstX = std::min(cX,cX+(xsize-1)*AXX);
    firstY = std::min(cY,cY+(ysize-1)*AYY);

    if (verbose>0)
      {
	double extendx = xsize*AXX;
	double extendy = ysize*AYY;
	double originx = std::min(oX,oX+extendx);
	double originy = std::min(oY,oY+extendy);
	std::cout << "GeoTIFF file " << filename << " opened" << std::endl;
	std::cout << " size=" << xsize << "x" << ysize << " with " << count << " band(s)" << std::endl;
	std::cout << " origin=" << originx << "," << originy << " extend=" << std::abs(extendx) << "," << std::abs(extendy)
		  << " first pixel at " << std::min(cX,cX+(xsize-1)*AXX) << "," << std::min(cY,cY+(ysize-1)*AYY) << std::endl;
	std::cout << "pixel width is " << std::setw(24)
				  << std::scientific
				  << std::showpoint
				  << std::setprecision(16) << AXX << " x " << AYY << std::defaultfloat << std::endl;
	if (std::abs(cX-round(cX))<1e-6 && std::abs(cY-round(cY))<1e-6)
	  std::cout << " vertex centered raster image" << std::endl;
	else
	  std::cout << " element centered raster image" << std::endl;	  
      }

    std::map<std::size_t,std::string> data_type_to_string;
    data_type_to_string[GDT_Unknown] = "unknown";
    data_type_to_string[GDT_Byte] = "Byte";
    data_type_to_string[GDT_UInt16] = "UInt16";
    data_type_to_string[GDT_Int16] = "Int16";
    data_type_to_string[GDT_UInt32] = "UInt32";
    data_type_to_string[GDT_Int32] = "Int32";
    data_type_to_string[GDT_Float32] = "Float32";
    data_type_to_string[GDT_Float64] = "Float64";
    data_type_to_string[GDT_CInt16] = "CInt16";
    data_type_to_string[GDT_CInt32] = "CInt32";
    data_type_to_string[GDT_CFloat32] = "CFloat32";
    data_type_to_string[GDT_CFloat64] = "CFloat64";
    
    // read raster band
    GDALRasterBand* band = poDataset->GetRasterBand(useband);
    dt = band->GetRasterDataType();
    _noDataValue = static_cast<T>(band->GetNoDataValue(&_hasNoDataValue));
    if (verbose>1)
      {
	std::cout << " reading band " << useband << ", DataType=" << dt << " (" << data_type_to_string[dt] << ")";
	if (_hasNoDataValue)
	  std::cout << " no data value=" << _noDataValue;
	else
	  std::cout << " no no data value!";
	std::cout << std::endl;
      }
    
    // now check if the data type in the image concides with the template parameter
    if ( band->GetRasterDataType()==GDT_Byte && std::is_same<T,GByte>::value==false)
      {
	GDALClose(poDataset);
	throw GeoTIFFReaderException("data type is not GByte");
      }
    if ( band->GetRasterDataType()==GDT_UInt16 && std::is_same<T,GUInt16>::value==false)
      {
	GDALClose(poDataset);
	throw GeoTIFFReaderException("data type is not GUInt16");
      }
    if ( band->GetRasterDataType()==GDT_Int16 && std::is_same<T,GInt16>::value==false)
      {
	GDALClose(poDataset);
	throw GeoTIFFReaderException("data type is not GInt16");
      }
    if ( band->GetRasterDataType()==GDT_UInt32 && std::is_same<T,GUInt32>::value==false)
      {
	GDALClose(poDataset);
	throw GeoTIFFReaderException("data type is not GUInt32");
      }
    if ( band->GetRasterDataType()==GDT_Int32 && std::is_same<T,GInt32>::value==false)
      {
	GDALClose(poDataset);
	throw GeoTIFFReaderException("data type is not GInt32");
      }
    if ( band->GetRasterDataType()==GDT_Float32 && std::is_same<T,float>::value==false)
      {
	GDALClose(poDataset);
	throw GeoTIFFReaderException("data type is not GFloat32");
      }
    if ( band->GetRasterDataType()==GDT_Float64 && std::is_same<T,double>::value==false)
      {
	GDALClose(poDataset);
	throw GeoTIFFReaderException("data type is not GFloat64");
      }

    // read block information
    int nXBlockSize, nYBlockSize;
    band->GetBlockSize( &nXBlockSize, &nYBlockSize );
    if (verbose>1)
      std::cout << " nXBlockSize=" << nXBlockSize << " nYBlockSize=" << nYBlockSize << std::endl;
    int nXBlocks = (band->GetXSize() + nXBlockSize - 1) / nXBlockSize;
    int nYBlocks = (band->GetYSize() + nYBlockSize - 1) / nYBlockSize;
    if (verbose>1)
      std::cout << " nXBlocks=" << nXBlocks << " nYBlocks=" << nYBlocks << std::endl;
    p = std::vector<T>(xsize*ysize); // allocate image
    T *pabyData = new T[nXBlockSize*nYBlockSize]; // allocate transfer buffer

    // loop over blocks and copy data
    for ( int iYBlock = 0; iYBlock < nYBlocks; iYBlock++ )
      for ( int iXBlock = 0; iXBlock < nXBlocks; iXBlock++ )
        {
          auto x=band->ReadBlock( iXBlock, iYBlock, pabyData );
          int nXValid, nYValid;
          band->GetActualBlockSize(iXBlock, iYBlock, &nXValid, &nYValid); // infer valid part of buffer
          for (int iy=0; iy<nYValid; iy++)
            { // loop over valid pixels in buffer
              auto iImage = ((iYBlock*nYBlockSize+iy)*xsize)+(iXBlock*nXBlockSize); // first pixel of line in image
              auto iBuffer = iy*nXBlockSize; // first pixel of line in buffer
              for (int ix=0; ix<nXValid; ix++) p[iImage++] = pabyData[iBuffer++];
            }
        }
    delete[] pabyData; // delete buffer

    // close data set
    GDALClose(poDataset);
  }

  GDALDataType getGDALDataType ()
  {
    return dt;
  }
  
  //! apply affine GeoTransform to arbitrary vector data type
  template<typename U>
  U transform (const U& x) const
  {
    U y(x); // make a copy
    y[0] =  oX + AXX*x[0] + AXY*x[1];
    y[1] =  oY + AYX*x[0] + AYY*x[1];
    return y;
  }
  
  //! Pixel size in longitude direction
  double originLong () const
  {
    return firstX;
  }

  //! Pixel size in latitude direction
  double originLat () const
  {
    return firstY;
  }

  //! number of pixels in longitude direction
  int sizeLong () const
  {
    return xsize;
  }

  //! number of pixels in latitude direction
  int sizeLat () const
  {
    return ysize;
  }

  //! longitude of pixel *center*
  double longitude (int x, int y) const
  {
    return cX + AXX*x + AXY*y;
  }

  //! latitude of pixel *center*
  double latitude (int x, int y) const
  {
    return cY + AYX*x + AYY*y;
  }

  //! Pixel size in longitude direction
  double dLong () const
  {
    return AXX;
  }

  //! Pixel size in latitude direction
  double dLat () const
  {
    return AYY;
  }

  // simplifed access operator with round brackets
  const T& operator() (int x, int y) const
  {
    return p[y*xsize+x];
  }
  T& operator() (int x, int y)
  {
    return p[y*xsize+x];
  }

  std::string name () const
  {
    return image_name;
  }
  
  // export underlying vector for reading
  const std::vector<T>& data ()
  {
    return p;
  }

  bool hasNoDataValue () const
  {
    return (_hasNoDataValue!=0);
  }

  T noDataValue () const
  {
    return _noDataValue;
  }
  
private:
  int xsize, ysize; // number of pixels in x and y
  int count; // number of raster bands in the data set
  double oX,oY;     // origin of image
  double cX,cY;     // origin of image offset by a half pixel
  double AXX, AXY, AYX, AYY; // transformation matrix
  std::string image_name; // name of the image
  std::vector<T> p; // image buffer
  double firstX, firstY; // position of first pixel in lower left corner (instead of up and left)!
  GDALDataType dt;
  int _hasNoDataValue;
  T _noDataValue;
};

#endif // Udune_hydro_HH
