#ifndef Udune_hydro_rasterdataset_HH
#define Udune_hydro_rasterdataset_HH

#include <exception>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

//! \brief Exception thrown by RasterDataSet
class RasterDataSetException : public std::exception
{
public:
  virtual const char* what() const throw()
  {
    return message.c_str();
  }

  RasterDataSetException (std::string m) : message(m) {}
  
private:
  std::string message;
};


  /*! \file
   *  \brief This file implements RasterDataSet, a georeferenced image
   */

  /** \brief Georeferenced image data
   *
   * \tparam T is the type used for each pixel
   */
template<typename T>
class RasterDataSet {
  double _originLong; // longitude of the lower left pixel in the image
  double _originLat;  // latitude of the lower left pixel in the image 
  double _dLong;      // width of a pixel in longitude direction
  double _dLat;       // width of a pixel in latitude direction
  int _sizeLong;      // number of pixels in longitude direction
  int _sizeLat;       // number of pixels in latitude direction
  T _fill_value;       // value used for initialization and missing pixels
  int verbose;        // verbose level, 0=silent, 1=overview, 2=detail
  std::vector<T> _data;// pixel data
  bool _hasFillValue;

public:
  // ==== expoorted types
  typedef T pixel_type; // export type for pixel
  typedef T value_type; // export type for pixel
  
  // ==== the constructors
  
  //! \brief Default constructor
  RasterDataSet () {}

  /**
   * \brief construct initialized RasterDataSet
   *
   * A RasterDataSet is an image with associated longitude/latitude coordinates for each pixel.
   * The coordinate of a pixel is the position of the center of the pixel on the globe.
   * RasterDataSets are always upright. The origin of the image is in the lower left corner.
   * Pixels are stored in lexicographic order from lower left to upper right.
   *
   * \param originLong_ longitude of the lower left pixel in the image
   * \param originLat_ latitude of the lower left pixel in the image
   * \param dLong_ width of a pixel in longitude direction
   * \param dLat_ width of a pixel in latitude direction
   * \param sizeLong_ number of pixels in longitude direction
   * \param sizeLat_ number of pixels in latitude direction
   * \param fill_value_ initial value for each pixel
   * \param verbose_ verbose level, 0=silent, 1=overview, 2=detail
   */
  RasterDataSet (double originLong_, double originLat_, double dLong_, double dLat_, int sizeLong_, int sizeLat_, T fill_value_, int verbose_=0)
    : _originLong(originLong_), _originLat(originLat_), _dLong(dLong_), _dLat(dLat_), _sizeLong(sizeLong_), _sizeLat(sizeLat_), _fill_value(fill_value_), verbose(verbose_),
      _data(sizeLong_*sizeLat_,fill_value_), _hasFillValue(true)
  {}

    /**
   * \brief construct initialized RasterDataSet
   *
   * A RasterDataSet is an image with associated longitude/latitude coordinates for each pixel.
   * The coordinate of a pixel is the position of the center of the pixel on the globe.
   * RasterDataSets are always upright. The origin of the image is in the lower left corner.
   * Pixels are stored in lexicographic order from lower left to upper right.
   *
   * \param originLong_ longitude of the lower left pixel in the image
   * \param originLat_ latitude of the lower left pixel in the image
   * \param dLong_ width of a pixel in longitude direction
   * \param dLat_ width of a pixel in latitude direction
   * \param sizeLong_ number of pixels in longitude direction
   * \param sizeLat_ number of pixels in latitude direction
   * \param t pass in a data array
   * \param verbose_ verbose level, 0=silent, 1=overview, 2=detail
   */
  RasterDataSet (double originLong_, double originLat_, double dLong_, double dLat_, int sizeLong_, int sizeLat_, const std::vector<T>& t, int verbose_=0)
    : _originLong(originLong_), _originLat(originLat_), _dLong(dLong_), _dLat(dLat_), _sizeLong(sizeLong_), _sizeLat(sizeLat_), verbose(verbose_),
      _data(t), _hasFillValue(false), _fill_value(0)
  {
    if (_data.size()!=_sizeLong*_sizeLat) _data.resize(_sizeLong*_sizeLat);
  }

  // ==== query methods

  /** \brief apply affine linear transformation from pixel space to long/lat coordinates
   * \tparam U arbitrary type with operator[] access. Same type is returned
   * \param x pixel position as arbitratry vector type.
   */
  template<typename U>
  U transform (const U& x) const
  {
    U y(x); // make a copy
    y[0] =  _originLong + _dLong*x[0];
    y[1] =  _originLat  + _dLat *x[1];
    return y;
  }
  
  //! \brief longitude of lower left pixel
  double originLong () const
  {
    return _originLong;
  }

  //! \brief latitude of lower left pixel
  double originLat () const
  {
    return _originLat;
  }
  
  //! \brief number of pixels in x (longitude)
  int sizeLong () const
  {
    return _sizeLong;
  }

  //! \brief number of pixels in y (latitude)
  int sizeLat () const
  {
    return _sizeLat;
  }

  //! \brief longitude of pixel center
  double longitude (int x, int y) const
  {
    return _originLong + _dLong*x;
  }

  //! \brief latitude of pixel center
  double latitude (int x, int y) const
  {
    return _originLat + _dLat*y;
  }

  //! \brief Pixel size in longitude
  double dLong () const
  {
    return _dLong;
  }

  //! \brief Pixel size in latitude
  double dLat () const
  {
    return _dLat;
  }

  // ==== access to pixel data

  //! \brief access with pixel coordinates
  const T& operator() (int x, int y) const
  {
    return _data[y*_sizeLong+x];
  }
  
  //! \brief const access with pixel coordinates
  T& operator() (int x, int y)
  {
    return _data[y*_sizeLong+x];
  }

  //! \brief const access with long/lat coordinates; this is read only and returns a value!
  T operator() (double lon, double lat) const
  {
    int x = round((lon-_originLong)/_dLong);
    int y = round((lat-_originLat)/_dLat);
    int i = y*_sizeLong+x;
    if (i>=0 && i<_data.size())
      return _data[i];
    else
      return _fill_value;
  }
  
  //! export underlying std::vector for reading
  const std::vector<T>& data () const
  {
    return _data;
  }

  //! export underlying std::vector for writing
  std::vector<T>& data ()
  {
    return _data;
  }

  bool hasNoDataValue () const
  {
    return _hasFillValue;
  }

  T noDataValue () const
  {
    return _fill_value;
  }

  // ==== higher level operations on the image

  /** \brief paste onother image into this image at the appropriate place
   *
   * Paste image into this image. Currently only images with the same pixel size
   * can be pasted. Computes the intersection of both images and pastes only the
   * relevant part. The input image may also have its origin in any corner. 
   *
   * \tparam Image type supporting the same operations as RasterDataSet; may also be a GeoTIFFImage
   *
   * \param image input image
   */
  template<typename Image>
  void paste (const Image& image)
  {
    // check if pixel width is the same
    if (std::abs(std::abs(image.dLong())-std::abs(_dLong))>1e-6 || std::abs(std::abs(image.dLat())-std::abs(_dLat))>1e-6)
      throw RasterDataSetException("pixel width does not match in paste, nothing done");

    // extent of image in longitude
    double minLong = image.longitude(0,0);
    double maxLong = image.longitude(image.sizeLong()-1,0);
    int flipLong = 1;
    if (minLong>maxLong) {
      std::swap(minLong,maxLong);
      flipLong = -1;
    }

    // extent of image in latitude
    double minLat = image.latitude(0,0);
    double maxLat = image.latitude(0,image.sizeLat()-1);
    int flipLat = 1;
    if (minLat>maxLat) {
      std::swap(minLat,maxLat);
      flipLat = -1;
    }

    // compute offset between image and canvas
    // offset_Long,offset_Lat is the offset of the origin of the image with respect to the origin of the canvas
    // i.e. origin_canvas + offset = origin_image (in pixel coordinates of the canvas)
    // ==> origin_image - offset = origin_canvas (in pixel coordinates of the image)
    int offset_Long = round((minLong - _originLong)/std::abs(image.dLong()));
    int offset_Lat = round((minLat - _originLat)/std::abs(image.dLat()));

    // compute the size of the common image area
    int sLong = std::min(offset_Long+image.sizeLong(),_sizeLong) - std::max(offset_Long,0);
    int sLat = std::min(offset_Lat+image.sizeLat(),_sizeLat) - std::max(offset_Lat,0);

    // compute origin of the copy area in the canvas
    int start_canvas_Long = std::max(offset_Long,0);
    int start_canvas_Lat = std::max(offset_Lat,0);

    // compute origin of the copy area in the image
    int start_image_Long = std::max(-offset_Long,0);
    int start_image_Lat = std::max(-offset_Lat,0);

    // now copy the intersection
    int OLong=0, OLat=0; // for flipping transformation
    if (flipLong<0) OLong = image.sizeLong()-1;
    if (flipLat<0) OLat = image.sizeLat()-1;
    if (verbose>0)
      {
        std::cout << "pasting image, copy area size " << sLong << "x" << sLat;
        std::cout << " from " << start_image_Long << "," << start_image_Lat << " to " << start_canvas_Long << "," << start_canvas_Lat;
        std::cout << " flip " << flipLong << "," << flipLat << std::endl;
      }
    if (image.hasNoDataValue())
      {
	auto no_data_value = image.noDataValue();
	if (sLong*sLat>0)
	  for (int j=0; j<sLat; ++j)
	    {
	      int cj = start_canvas_Lat+j;
	      int ij = OLat+flipLat*(start_image_Lat+j);
	      for (int i=0; i<sLong; ++i)
		{
		  if (image(OLong+flipLong*(start_image_Long+i),ij)==no_data_value)
		    _data[cj*_sizeLong+start_canvas_Long+i] = _fill_value;
		  else
		    _data[cj*_sizeLong+start_canvas_Long+i] =  image(OLong+flipLong*(start_image_Long+i),ij);
		  if (std::is_same<T,float>::value==true || std::is_same<T,double>::value==true)
		    if (isnan(_data[cj*_sizeLong+start_canvas_Long+i]))
		      _data[cj*_sizeLong+start_canvas_Long+i] = _fill_value;
		}
	    }
      }
    else
      {
	if (sLong*sLat>0)
	  for (int j=0; j<sLat; ++j)
	    {
	      int cj = start_canvas_Lat+j;
	      int ij = OLat+flipLat*(start_image_Lat+j);
	      for (int i=0; i<sLong; ++i)
		_data[cj*_sizeLong+start_canvas_Long+i] =  image(OLong+flipLong*(start_image_Long+i),ij);
	    }
      }      
  }

};

#endif
