// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <cmath>
#include <string>  
#include <iostream> 
#include <sstream>
// C++ includes
#include<math.h>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif

//GDAL includes. We require that GDAL is found
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()


/** list all tiles needed for a certain region
    The arguments are given in latitude and longitude:
    x_lower   longitude of most west point, west of Greenwich is negative, east is positive
    y_lower   latitude of most south point, south of equator is negative, north of equator is positive
    x_upper   longitude of most east point, west of Greenwich is negative, east is positive
    y_upper   latitude of most north point, south of equator is negative, north of equator is positive
 */
void listTiles (double x_lower, double y_lower, double x_upper, double y_upper)
{
  // tiles have the format [n|s]YY_[e,w]XXX_1arc_v3.tif
  int xint_lower = floor(x_lower);
  int yint_lower = floor(y_lower);
  int xint_upper = ceil(x_upper);
  int yint_upper = ceil(y_upper);
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

/** check if all tiles are available
 */
bool checkTiles ( double x_lower, double y_lower, double x_upper, double y_upper, std::string path)
{
}

/** A class that reads a GeoTIFFImage
    
    - holds the data in a buffer and lets you access the pixels
    - provides size of the image
    - provides positions of pixel centers in lat-long
    - Note: images are y-downwards!
 */
class GeoTIFFImage {
public:
  GeoTIFFImage (std::string filename, std::string path="")
  {
    GDALDataset  *poDataset;
    poDataset = (GDALDataset *) GDALOpen( filename.c_str(), GA_ReadOnly );
    if ( poDataset == NULL )
      {
        std::cout << "GDALOpen() failed" << std::endl;
        // Need to throw an exception here
      }

    // read size
    xsize = poDataset->GetRasterXSize();
    ysize = poDataset->GetRasterYSize();

    // read transformation to lat long coordinates
    double GT[6];
    if ( poDataset->GetGeoTransform( GT ) != CE_None )
      {
        std::cout << "GetGeoTransform() failed" << std::endl;
        // Need to throw an exception here
      }
    oX = GT[0];  oY = GT[3];
    AXX = GT[1]; AXY = GT[2];
    AYX = GT[4]; AYY = GT[5];
    
    // read file contents
    GDALRasterBand* band = poDataset->GetRasterBand(1);
    int nXBlockSize, nYBlockSize;
    band->GetBlockSize( &nXBlockSize, &nYBlockSize );
    int nXBlocks = (band->GetXSize() + nXBlockSize - 1) / nXBlockSize;
    int nYBlocks = (band->GetYSize() + nYBlockSize - 1) / nYBlockSize;
    p = new GInt16[xsize*ysize]; // allocate image
    GInt16 *pabyData = new GInt16[nXBlocks*nYBlocks]; // allocate transfer buffer

    // loop over blocks
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
  }

  GeoTIFFImage (const GeoTIFFImage& image)
  {
    // deep copy
    xsize = image.xsize;
    ysize = image.ysize;
    
    oX = image.oX;
    oY = image.oY;
    AXX = image.AXX;
    AXY = image.AXY;
    AYX = image.AYX;
    AYY = image.AYY;
    
    p = new GInt16[xsize*ysize]; // allocate image
    for (int i=0; i<xsize*ysize; i++) p[i] = image.p[i];
  }
  
  // destructor
  ~GeoTIFFImage ()
  {
    delete[] p;
  }
  
  // number of pixels in x (longitude)
  int sizeX () const
  {
    return xsize;
  }

  // number of pixels in y (latitude)
  int sizeY () const
  {
    return ysize;
  }

  // longitude of pixel center
  double longitude (int x, int y) const
  {
    return oX + AXX*(x+0.5) + AXY*(y+0.5);
  }

  // latitude of pixel center
  double latitude (int x, int y) const
  {
    return oY + AYX*(x+0.5) + AYY*(y+0.5);
  }

  // read pixel 
  GInt16 getPixel (int x, int y) const
  {
    return p[y*xsize+x];
  }

  // read pixel 
  void putPixel (int x, int y, GInt16 value)
  {
    p[y*xsize+x] = value;
  }

  // simplifed access operator with round brackets
  const GInt16& operator() (int x, int y) const
  {
    return p[y*xsize+x];
  }
  GInt16& operator() (int x, int y)
  {
    return p[y*xsize+x];
  }

private:
  int xsize, ysize; // number of pixels in x and y
  double oX,oY;     // origin of image
  double AXX, AXY, AYX, AYY; // transformation matrix
  GInt16 *p; // image buffer
};

int driver (int argc, char** argv)
{
  GDALDataset  *poDataset;
  GDALAllRegister();
  std::string name(argv[1]); // "n49_e009_1arc_v3.tif"
  poDataset = (GDALDataset *) GDALOpen( name.c_str(), GA_ReadOnly );
  if( poDataset == NULL )
    {
      std::cout << "GDALOpen() failed" << std::endl;
      return 1;
    }
  
  // getting meta data
  std::cout << poDataset->GetDriver()->GetDescription()
            << poDataset->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME )
            << std::endl;
  auto xsize = poDataset->GetRasterXSize();
  auto ysize = poDataset->GetRasterYSize();
  std::cout << "xsize=" << xsize
            << " ysize=" << ysize
            << " bands=" << poDataset->GetRasterCount()
            << std::endl;
  GDALRasterBand* band = poDataset->GetRasterBand(1);
  std::cout << "data type is " << GDALGetDataTypeName(band->GetRasterDataType()) << std::endl;
  if( poDataset->GetProjectionRef()  != NULL )
    std::cout << "Projection is " << poDataset->GetProjectionRef() << std::endl;
  double GT[6];
  
  if( poDataset->GetGeoTransform( GT ) != CE_None )
    {
      std::cout << "GetGeoTransform() failed" << std::endl;
      return 1;
    }
  // There are 3601x3601 pixels
  double DX=GT[1], DY=GT[5]; // extend of a pixel in x (longitued) and y (latitude)
  double X0=GT[0], Y0=GT[3]; // position of origin of the tile which is at the corner of a pixel
  std::cout << "Pixel Size = " << DX << "," << DY << std::endl;
  std::cout << "Origin = " << X0 << "," << Y0 << std::endl;
  // The center of the upper left pixel is 
  std::cout << "(0,0) -> " << X0 + 0.5*DX << "," << Y0 + 0.5*DY << std::endl;
  // The center of the lower right pixel is 
  std::cout << "(" << xsize-1 << ","
            << ysize-1 << ") -> "
            << X0 + (xsize-0.5)*DX  << ","
            << Y0 + (ysize-0.5)*DY
            << std::endl;
  // this means the tiles are actually overlapping by one pixel, i.e. in each row,
  // the first pixel is the last pixel from the previous tile and the same for the columns 

  /* Transformation:
     1) Pixel data (i,j), one tile is (0,0) to (3600,3600)
     These are actually the the values in the centers of the pixels. And the tiles are overlapping
     2) Latidute, longitude (phi,lambda) are provided by the affine geometry transformation
     3a) Position in meters in a flat model
     3b) Position in a whole earth model 
     WGS84 spheroid at sea level:
     one latitudinal second measures 30.715 meters, 
     one latitudinal minute is 1843  meters and 
     one latitudinal degree is 110.6 kilometres
     On the WGS84 spheroid, the length in meters of a degree of latitude at latitude φ (that is, the number of meters you would have to travel along a north–south line to move 1 degree in latitude, when at latitude φ radians), is about

     111132.92 − 559.82 cos ⁡ 2 φ + 1.175 cos ⁡ 4 φ − 0.0023 cos ⁡ 6 φ {\displaystyle 111132.92-559.82\,\cos 2\varphi +1.175\,\cos 4\varphi -0.0023\,\cos 6\varphi } 111132.92-559.82\,\cos 2\varphi +1.175\,\cos 4\varphi -0.0023\,\cos 6\varphi [10]

     The returned measure of meters per degree latitude varies continuously with latitude.

     Similarly, the length in meters of a degree of longitude can be calculated as

     111412.84 cos ⁡ φ − 93.5 cos ⁡ 3 φ + 0.118 cos ⁡ 5 φ {\displaystyle 111412.84\,\cos \varphi -93.5\,\cos 3\varphi +0.118\,\cos 5\varphi } {\displaystyle 111412.84\,\cos \varphi -93.5\,\cos 3\varphi +0.118\,\cos 5\varphi }[10]
  */

    
  /* next steps:
     - make a seperate application, not so important
     - connect data with a yasp grid of appropriate size (in meters) 
     - output to paraview
     - think about an interface for data reading: lat,lat long region
     - data 
  */
}

int main(int argc, char** argv)
{
  try{
    // Maybe initialize MPI
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    std::cout << "Hello World! This is dune-hydro." << std::endl;
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()
               <<" processes!"<<std::endl;

    // register GDAL drivers once at the beginning
    GDALAllRegister();

    // test file lister
    listTiles(7.7, 48.0, 8.5, 50.0);

    // load an image
    GeoTIFFImage image("n49_e008_1arc_v3.tif");

    // make a yaspgrid
    const int dim = 2;
    using Grid = Dune::YaspGrid<dim>; 
    Dune::FieldVector<double,dim> len;
    len[0]=std::abs(image.longitude(image.sizeX()-1,0)-image.longitude(0,0))*3600*30.715; 
    len[1]=std::abs(image.latitude(0,image.sizeY()-1)-image.latitude(0,0))*3600*30.715;
    std::cout << "domain extend:" << len[0] << " " << len[1] << std::endl;
    std::array<int,dim> cells;
    cells[0]=image.sizeX()-1;
    cells[1]=image.sizeY()-1;
    if (cells[0]!=3600 || cells[1]!=3600)
      {
        std::cout << "only works with 3600^2 cells" << std::endl;
        return 1;
      }
    cells[0] /= 16; cells[1] /= 16; 
    std::cout << "cells:" << cells[0] << " " << cells[1] << std::endl;
    Grid grid(len,cells);
    grid.globalRefine(4);

    // allocate a hierarchy of representations
    int L = grid.maxLevel();
    std::vector<std::vector<float>> heightdata(L+1);
    for (int l=L; l>=0; l--) heightdata[l].resize(grid.levelGridView(l).indexSet().size(0));
    std::vector<std::vector<float>> curvature(L+1);
    for (int l=L; l>=0; l--) curvature[l].resize(grid.levelGridView(l).indexSet().size(0));
    std::vector<std::vector<float>> gradientn(L+1);
    for (int l=L; l>=0; l--) gradientn[l].resize(grid.levelGridView(l).indexSet().size(0));

    // fill data on finest level
    {
      auto gv = grid.leafGridView();
      auto& indexset = gv.indexSet();
      int xsize=image.sizeX()-1;
      int ysize=image.sizeY()-1;
      for (int y=0; y<ysize; y++)
        for (int x=0; x<xsize; x++)
          heightdata[L][y*xsize+x] = image(x,xsize-1-y);
      Dune::VTKWriter<decltype(gv)> vtkwriter(gv,Dune::VTK::conforming);
      vtkwriter.addCellData(heightdata[L],"height");
      vtkwriter.write("leafvtk",Dune::VTK::appendedraw);
    }
    
    // compute coarser representations of height
    for (int l=L; l>0; l--)
      {
        // level wise access to elements
        auto finegv = grid.levelGridView(l);
        auto coarsegv = grid.levelGridView(l-1);
        auto& fineindexset = finegv.indexSet();
        auto& coarseindexset = coarsegv.indexSet();
        for (int i=0; i<coarseindexset.size(0); i++)
          heightdata[l-1][i] = 0.0; // initialize coarse grid
        for (const auto& e : elements(finegv))
          heightdata[l-1][coarseindexset.index(e.father())] += 0.25*heightdata[l][fineindexset.index(e)];
      }

    // compute curvature on each level
    for (int l=L; l>=0; l--)
      {
        float s = pow(2.0,L)*pow(2.0,L);
        auto gv = grid.levelGridView(l);
        auto& indexset = gv.indexSet();
        for (const auto& e : elements(gv)) {
          float c = 0.0;
          for (const auto& is : intersections(gv,e))
            if (is.neighbor())
              c += heightdata[l][indexset.index(is.outside())]-heightdata[l][indexset.index(e)];
          curvature[l][indexset.index(e)] = s*c;
        }
      }

    // compute gradient on each level
    for (int l=L; l>=0; l--)
      {
        float s = pow(2.0,L);
        auto gv = grid.levelGridView(l);
        auto& indexset = gv.indexSet();
        for (const auto& e : elements(gv)) {
          float c = heightdata[l][indexset.index(e)];
          float n[4]; int i=0;
          for (const auto& is : intersections(gv,e)) {
            n[i] = c;
            if (is.neighbor()) n[i] = heightdata[l][indexset.index(is.outside())];
            i++;
          }
          gradientn[l][indexset.index(e)] = s*(std::abs(n[1]-c)+std::abs(n[3]-c));
        }
      }

    // write vtk output
    for (int l=0; l<=L; l++)
      {
        auto finegv = grid.levelGridView(l);
        Dune::VTKWriter<decltype(finegv)> vtkwriter(finegv,Dune::VTK::conforming);
        vtkwriter.addCellData(heightdata[l],"height");
        vtkwriter.addCellData(curvature[l],"curvature");
        vtkwriter.addCellData(gradientn[l],"gradientn");
        std::stringstream fname;
        fname << "vtkfile" << l;      
        vtkwriter.write(fname.str(),Dune::VTK::appendedraw);
      }

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
