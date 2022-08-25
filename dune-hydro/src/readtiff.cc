// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/*************************************/
/* This example solves surface flow in the vicinity of Heidelberg
   with 30m resolution and constant rainfall. Uses structured mesh
   and finite volume method.
 */
/*************************************/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <cmath>
#include <string>  
#include <iostream> 
#include <sstream>
#include <queue>
#include <algorithm>
// C++ includes
#include<math.h>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include<dune/common/exceptions.hh> // We use exceptions
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/io/file/vtk.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
// dune-pdelab includes
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/variablefactories.hh>
#include<dune/pdelab/localoperator/l2.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/newton/newton.hh>
#include<dune/pdelab/solver/newton.hh>

// dune-vtk includes
#include <dune/vtk/writers/vtkimagedatawriter.hh>
#include <dune/vtk/writers/vtkrectilineargridwriter.hh>
#include <dune/vtk/writers/vtkstructuredgridwriter.hh>
#include <dune/vtk/writers/vtkunstructuredgridwriter.hh>
#include <dune/vtk/datacollectors/yaspdatacollector.hh>
#include <dune/vtk/datacollectors/spdatacollector.hh>

// include stuff from dune-hydro
#include<dune/hydro/nonlineardiffusionfv.hh>
#include<dune/hydro/nonlineardiffusionfv_velocity.hh>
#include<dune/hydro/geotiffreader.hh>
#include<dune/hydro/netcdfreader.hh>
#include<dune/hydro/rasterdataset.hh>


// classification of segments
enum SegmentType { unknown, basin, plain };

// data to be stored for each link from one segment to another
struct LinkData
{
  int passheight; // height of lowest pixel between two segments
  int index; // index of pixel at the passheight within the segment
  LinkData () : passheight(1000000) {}
  LinkData (int ph, int i)
    : passheight(ph), index(i)
  {}
};

// data to be stored for each segment
struct SegmentData
{
  int seedindex;
  SegmentType type;
  int size;
  int drainto;
  int accumulation;
  double average_height;
  bool boundary;
  std::map<int,LinkData> links;
  SegmentData ()
    : seedindex(-1), type(unknown), size(0), drainto(-1), accumulation(0), boundary(false)
  {}
  SegmentData (int i, SegmentType t, int s)
    : seedindex(i), type(t), size(s), drainto(-1), accumulation(0), boundary(false)
  {}
};

template<typename T>
RasterDataSet<T> arithmetic_average (RasterDataSet<T>& image_in, int diameter)
{
  // prepare result
  RasterDataSet<T> image_out(image_in);
  std::vector<T>& out = image_out.data();

  // analyse size etc of input image
  auto& in=image_in.data(); // access to raw image data
  int N=in.size(); // number of pixels in image
  int m=image_in.sizeLong(); // pixels per row
  int n=image_in.sizeLat(); // pixels per column

  // compute average value of pixels
  for (int j=0; j<n; j++)
    for (int i=0; i<m; i++)
      {
        double sum=0.0;
        for (int jj=std::max(j-diameter,0); jj<std::min(j+diameter+1,n); jj++)
          for (int ii=std::max(i-diameter,0); ii<std::min(i+diameter+1,m); ii++)
            sum += ((double)in[jj*m+ii]);
        out[j*m+i] = sum/((std::min(j+diameter+1,n)-std::max(j-diameter,0))*(std::min(i+diameter+1,m)-std::max(i-diameter,0)));
      }
  return image_out;
}

template<typename T>
RasterDataSet<T> geometric_average (RasterDataSet<T>& image_in, int diameter)
{
  // prepare result
  RasterDataSet<T> image_out(image_in);
  std::vector<T>& out = image_out.data();

  // analyse size etc of input image
  auto& in=image_in.data(); // access to raw image data
  int N=in.size(); // number of pixels in image
  int m=image_in.sizeLong(); // pixels per row
  int n=image_in.sizeLat(); // pixels per column

  // compute average value of pixels
  for (int j=0; j<n; j++)
    for (int i=0; i<m; i++)
      {
        double product=1.0;
        for (int jj=std::max(j-diameter,0); jj<std::min(j+diameter+1,n); jj++)
          for (int ii=std::max(i-diameter,0); ii<std::min(i+diameter+1,m); ii++)
            product *= ((double)in[jj*m+ii]);
        int npixel = (std::min(j+diameter+1,n)-std::max(j-diameter,0))*(std::min(i+diameter+1,m)-std::max(i-diameter,0));
        out[j*m+i] = pow(product,1.0/((double)npixel));
      }
  return image_out;
}

template<typename T>
RasterDataSet<T> median (RasterDataSet<T>& image_in, int diameter)
{
  // prepare result
  RasterDataSet<T> image_out(image_in);
  std::vector<T>& out = image_out.data();
  diameter = std::max(1,diameter);

  // analyse size etc of input image
  auto& in=image_in.data(); // access to raw image data
  int N=in.size(); // number of pixels in image
  int m=image_in.sizeLong(); // pixels per row
  int n=image_in.sizeLat(); // pixels per column

  // compute average value of pixels
  for (int j=0; j<n; j++)
    for (int i=0; i<m; i++)
      {
        std::vector<T> values;
        for (int jj=std::max(j-diameter,0); jj<std::min(j+diameter+1,n); jj++)
          for (int ii=std::max(i-diameter,0); ii<std::min(i+diameter+1,m); ii++)
            values.push_back(in[jj*m+ii]);
        std::sort(values.begin(),values.end());
        if (values.size()%2==0)
          out[j*m+i] = (values[values.size()/2]+values[values.size()/2-1])/2;
        else
          out[j*m+i] = values[values.size()/2];
      }
  return image_out;
}

template<typename T>
RasterDataSet<T> outlier (RasterDataSet<T>& image_in, int diameter, T threshold)
{
  // prepare result
  RasterDataSet<T> image_out(image_in);
  std::vector<T>& out = image_out.data();
  diameter = std::max(1,diameter);

  // analyse size etc of input image
  auto& in=image_in.data(); // access to raw image data
  int N=in.size(); // number of pixels in image
  int m=image_in.sizeLong(); // pixels per row
  int n=image_in.sizeLat(); // pixels per column

  // compute average value of pixels
  for (int j=0; j<n; j++)
    for (int i=0; i<m; i++)
      {
        std::vector<T> values;
        for (int jj=std::max(j-diameter,0); jj<std::min(j+diameter+1,n); jj++)
          for (int ii=std::max(i-diameter,0); ii<std::min(i+diameter+1,m); ii++)
            values.push_back(in[jj*m+ii]);
        std::sort(values.begin(),values.end());
        T median;
        if (values.size()%2==0)
          median = (values[values.size()/2]+values[values.size()/2-1])/2;
        else
          median = values[values.size()/2];
        if (std::abs(median-in[j*m+i])>=threshold)
          {
            std::cout << "i=" << i << " j=" << j << " value=" << in[j*m+i] << " median=" << median << std::endl;
            out[j*m+i] = median;
          }
        else
          out[j*m+i] = in[j*m+i];
      }
  return image_out;
}

template<typename T>
RasterDataSet<T> variance (RasterDataSet<T>& image_in, int diameter)
{
  // ensure pixel type is int
  if (std::is_same<T,int>::value==false)
    {
      std::cout << "pixel type must be int" << std::endl;
      exit(1);
    }

  // prepare result
  RasterDataSet<T> image_out(image_in);
  std::vector<T>& out = image_out.data();

  // analyse size etc of input image
  auto& in=image_in.data(); // access to raw image data
  int N=in.size(); // number of pixels in image
  int m=image_in.sizeLong(); // pixels per row
  int n=image_in.sizeLat(); // pixels per column

  // compute average value of pixels
  for (int j=0; j<n; j++)
    for (int i=0; i<m; i++)
      {
        double avg=0.0;
        for (int jj=std::max(j-diameter,0); jj<std::min(j+diameter+1,n); jj++)
          for (int ii=std::max(i-diameter,0); ii<std::min(i+diameter+1,m); ii++)
            avg += ((double)in[jj*m+ii]);
        int npixel = (std::min(j+diameter+1,n)-std::max(j-diameter,0))*(std::min(i+diameter+1,m)-std::max(i-diameter,0));
        avg /= npixel;
        double var=0.0;
        for (int jj=std::max(j-diameter,0); jj<std::min(j+diameter+1,n); jj++)
          for (int ii=std::max(i-diameter,0); ii<std::min(i+diameter+1,m); ii++)
            var += (((double)in[jj*m+ii])-avg)*(((double)in[jj*m+ii])-avg);
        var /= npixel;
        out[j*m+i] = var;
      }
  return image_out;
}


template<typename T>
RasterDataSet<T> flow_accumulation (RasterDataSet<T>& image)
{
  // ensure pixel type is int
  if (std::is_same<T,int>::value==false)
    {
      std::cout << "pixel type must be int" << std::endl;
      exit(1);
    }
    
  // analyse size etc of input image
  auto& data=image.data(); // access to raw image data
  int N=data.size(); // number of pixels in image
  int m=image.sizeLong(); // pixels per row
  int n=image.sizeLat(); // pixels per column
  std::cout << "N=" << N << " m=" << m << " n=" << n << std::endl;

  /* prepare cases
   *  7 6 8
   *  1 0 2
   *  4 3 5
   */
  std::vector<std::vector<int>> neighbors ={{-m-1,-m,-m+1,-1,1,m-1,m,m+1}, // 0
                                            {-m,-m+1,1,m,m+1}, // 1
                                            {-m-1,-m,-1,m-1,m}, // 2
                                            {-1,1,m-1,m,m+1}, // 3
                                            {1,m,m+1}, // 4
                                            {-1,m-1,m}, // 5
                                            {-m-1,-m,-m+1,-1,1}, // 6
                                            {-m,-m+1,1}, // 7
                                            {-m-1,-m,-1}}; // 8
  std::vector<std::vector<int>> neighbors5 ={{-m,-1,1,m}, // 0
                                             {-m,1,m}, // 1
                                             {-m,-1,m}, // 2
                                             {-1,1,m}, // 3
                                             {1,m}, // 4
                                             {-1,m}, // 5
                                             {-m,-1,1}, // 6
                                             {-m,1}, // 7
                                             {-m,-1}}; // 8

  // allocate output image 
  RasterDataSet<T> result_image(image);
  std::vector<T>& result = result_image.data();
  
  // compute accumulation and drainto fields
  std::vector<int> drainto(N,-1); // a pixel drains to the lowest neighbor that is strictly lower than itself; -1 means there is none
  std::vector<int> accumulation(N,1); // number of other pixels draining to the given pixel + 1 (for itself)
  // sort pixels by height
  std::cout << "sort" << std::endl;
  using X = std::pair<int,int>;
  std::vector<X> x; // pixels sorted by height 
  for (int j=0; j<n; j++)
    for (int i=0; i<m; i++)
      {
        int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0));
        int index = j*m+i;
        T value = data[index];
        x.push_back(X(value,index));
      }
  std::sort(x.begin(),x.end());
  // accumulate to lowest neighbor
  std::cout << "accumulate" << std::endl;
  for (auto it=x.rbegin(); it!=x.rend(); ++it) // this is important that we proceed from high to low
    {
      int height = it->first;
      int index = it->second;
      int i=index%m;
      int j=index/m;
      int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0));
      int minheight=height;
      int mink=-1;
      for (int k=0; k<neighbors[c].size(); k++)
        {
          int index_neighbor=index+neighbors[c][k];
          if (data[index_neighbor]<minheight)
            {
              minheight = data[index_neighbor];
              mink = k;
            }
        }
      if (mink>=0)
        {
          // we have a neighbor with smaller height; accumulate this pixel to the neighbor
          int index_neighbor=index+neighbors[c][mink];
          accumulation[index_neighbor] += accumulation[index];
          drainto[index] = index_neighbor;
        }
    }

  /*
   * pixel i is a pit, when (drainto[i]==-1 && accumulation[i]>1)
   * pixel i is plain, when (drainto[i]==-1 && accumulation[i]==1)
   */
  // first statistics
  int count = 0;
  for (int i=0; i<N; i++)
    if (drainto[i]<0 && accumulation[i]>1)
      count++;
  std::cout << count << " pits or out of " << N << " pixels (" << (int)round(100.0*count/N) << "%)" << std::endl;

  // create segments
  std::vector<int> segment(N,-1); // segment number for each pixel
  std::vector<SegmentData> segments;
  std::vector<int> distance(N,-1); // distance to the seed along the accumulation tree; flat segments are omitted first
  int assigned_pixels=0; // count pixels assigned to a segment
  
  // first round; classify the real pits
  for (auto it=x.begin(); it!=x.end(); ++it)
    if (drainto[it->second]==-1 && accumulation[it->second]>1)
    {
      int seedindex = it->second;
      // now we have no drainage and accumulate more than one pixel, i.e. a pit
      if (segment[seedindex]!=-1) std::cout << "Uuups : pit pixel should not be assigned yet" << std::endl;

      // start a new segment
      int segment_number = segments.size();
      segments.push_back(SegmentData(seedindex,basin,0));
      segments[segment_number].average_height = 0.0;
      std::queue<int> q;
      q.push(seedindex);
      while (!q.empty())
        {
          // get next pixel from queue
          int index = q.front(); q.pop();
          // if pixel is already visited then skip it
          if (segment[index]!=-1) continue;
          // otherwise assign it to the current segment
          segment[index] = segment_number;
          segments[segment_number].average_height = (segments[segment_number].average_height*segments[segment_number].size + data[index])/(1+segments[segment_number].size);
          segments[segment_number].size += 1;
          if (drainto[index]<0)
            distance[index] = 0; // this is the pit.
          else
            distance[index] = distance[drainto[index]]+1; // increase distance by one
          assigned_pixels++;
          // now go to neighbors and see if they are draining to this pixel
          int i=index%m;
          int j=index/m;
          int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0));
          if (c!=0) segments[segment_number].boundary = true;
          for (int k=0; k<neighbors[c].size(); k++)
            {
              int index_neighbor=index+neighbors[c][k];
              if (drainto[index_neighbor]==index)
                q.push(index_neighbor);
            }
        }
      // a segment has been finished
    }
  std::cout << segments.size() << " segments classified as basin, "
            << assigned_pixels << " out of " << N << " pixels assigned (" <<  (int)round(100.0*assigned_pixels/N) << "%)" << std::endl;

  // second round; classify the exactly flat plains (connected pixels of same height)
  int s1 = segments.size();
  for (auto it=x.begin(); it!=x.end(); ++it)
    if (drainto[it->second]==-1 && accumulation[it->second]==1 && segment[it->second]==-1)
    {
      // now we have a plain pixel not assigned yet
      int seedindex = it->second;

      // start a new segment
      int segment_number = segments.size();
      segments.push_back(SegmentData(seedindex,plain,0));
      segments[segment_number].average_height = data[seedindex];
      std::queue<int> q;
      q.push(seedindex);
      while (!q.empty())
        {
          // get next pixel from queue
          int index = q.front(); q.pop();
          // if pixel is already visited then skip it
          if (segment[index]!=-1) continue;
          // otherwise assign it to the current segment
          segment[index] = segment_number;
          segments[segment_number].size += 1;
          assigned_pixels++;
          // now go to neighbors and see if they are draining to this pixel
          int i=index%m;
          int j=index/m;
          int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0));
          if (c!=0) segments[segment_number].boundary = true;
          for (int k=0; k<neighbors[c].size(); k++)
            {
              int index_neighbor=index+neighbors[c][k];
              if (data[index_neighbor]==data[index])
                q.push(index_neighbor);
            }
        }
      // a segment has been finished
    }
  std::cout << segments.size()-s1 << " segments classified as plain, "
            << assigned_pixels << " out of " << N << " pixels assigned (" <<  (int)round(100.0*assigned_pixels/N) << "%)" << std::endl;
  if (assigned_pixels!=N)
    std::cout << "Error: all pixels should now be assigned" << std::endl;

  // build segment graph with pass heights
  std::cout << "build segment graph" << std::endl;
  for (int j=0; j<n; j++)
    for (int i=0; i<m; i++)
      {
        int index = j*m+i;
        int my_segment_number = segment[index];
        int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0));
        for (int k=0; k<neighbors[c].size(); k++)
          {
            int index_neighbor=index+neighbors[c][k];
            int neighbor_segment_number = segment[index_neighbor];
            if (my_segment_number!=neighbor_segment_number)
              {
                auto& mylinks = segments[my_segment_number].links;
                auto& nblinks = segments[neighbor_segment_number].links;
                int passheight = std::max(data[index],data[index_neighbor]);
                auto it = mylinks.find(neighbor_segment_number);
                if (it==mylinks.end())
                  {
                    // new neighborship that has not been discovered before
                    mylinks[neighbor_segment_number] = LinkData(passheight,index);
                    nblinks[my_segment_number] = LinkData(passheight,index_neighbor);
                  }
                else
                  // the two segments are neighbors already
                  if (passheight<(it->second).passheight)
                    {
                      (it->second).passheight = passheight;
                      (it->second).index = index;
                      nblinks[my_segment_number].passheight = passheight;
                      nblinks[my_segment_number].index = index_neighbor;
                    }
              }
          }
      }

  // // report segments
  // for (int i=0; i<segments.size(); ++i)
  //   xs.push_back(X(segments[i].size,i));
  // std::sort(xs.begin(),xs.end());
  // std::cout << " largest segments:" << std::endl;
  // int cnt=0;
  // for (auto it=xs.rbegin(); it!=x.rend(); ++it)
  //   {
  //     if (cnt>1) break;
  //     cnt++;
  //     auto& s = segments[it->second];
  //     if (s.type==basin)
  //       std::cout << " basin: number=" << it->second << " size=" << s.size << " seed=" << s.seedindex << " lowest=" << data[s.seedindex]
  //                 << " aheight=" << s.average_height << " bnd=" << s.boundary << std::endl;
  //     if (s.type==plain)
  //       std::cout << " plain: number=" << it->second << " size=" << s.size << " seed=" << s.seedindex
  //                 << " aheight=" << s.average_height << " bnd=" << s.boundary << std::endl;
  //     for (auto it = s.links.begin(); it!=s.links.end(); ++it)
  //       {
  //         std::cout << "   nbsegment=" << it->first << " passheight=" << it->second.passheight
  //                   << " lowest in nb=" << data[segments[it->first].seedindex]
  //                   << " size=" << segments[it->first].size
  //                   << " type=" << segments[it->first].type
  //                   << std::endl;
  //       }
  //   }
  // std::cout << " smallest segments:" << std::endl;
  // cnt=0;
  // for (auto it=xs.begin(); it!=x.end(); ++it)
  //   {
  //     if (cnt>20) break;
  //     cnt++;
  //     auto& s = segments[it->second];
  //     if (s.type==basin)
  //       std::cout << " basin: number=" << it->second << " size=" << s.size << " seed=" << s.seedindex << " lowest=" << data[s.seedindex]
  //                 << " aheight=" << s.average_height << " bnd=" << s.boundary << std::endl;
  //     if (s.type==plain)
  //       std::cout << " plain: number=" << it->second << " size=" << s.size << " seed=" << s.seedindex
  //                 << " aheight=" << s.average_height << " bnd=" << s.boundary << std::endl;
  //     for (auto it = s.links.begin(); it!=s.links.end(); ++it)
  //       {
  //         std::cout << "   nbsegment=" << it->first << " passheight=" << it->second.passheight
  //                   << " lowest in nb=" << data[segments[it->first].seedindex]
  //                   << " size=" << segments[it->first].size
  //                   << " type=" << segments[it->first].type
  //                   << std::endl;
  //       }
  //   }

  // draining of segments: low hanging fruits
  // sort segments by height of the seedindex
  std::vector<X> xs; // segments sorted by size
  for (int i=0; i<segments.size(); ++i)
    xs.push_back(X(data[segments[i].seedindex],i));
  std::sort(xs.begin(),xs.end());
  
  std::cout << "low hanging fruit: drain segments" << std::endl;
  int counter1=0;
  int counter2=0;
  for (auto it2=xs.rbegin(); it2!=xs.rend(); ++it2)
    {
      int mylow = it2->first;
      int mysegmentnumber = it2->second;
      auto& mysegment = segments[mysegmentnumber]; // this is our segment

      if (mysegment.links.size()==0)
        continue;

      if (mysegment.drainto!=-1)
        {
          std::cout << "segment drains already; this should not happen" << std::endl;
          continue;
        }

      // sort links by 1) passheight, 2) height of seedindex of neighbor
      using Y = std::pair<int,X>;
      std::vector<Y> sorted_links;
      for (auto lit=mysegment.links.begin(); lit!=mysegment.links.end(); ++lit)
        sorted_links.push_back(Y(lit->second.passheight,X(data[segments[lit->first].seedindex],lit->first)));
      std::sort(sorted_links.begin(),sorted_links.end());

      // criterion 1: no pass necessary
      // - pass height to neighbor is equal to our lowest point
      // - the lowest point in the neighbor is lower than our lowest point
      // - the neighbor segment does not drain to me
      for (auto& l : sorted_links)
        {
          int passheight = l.first;
          int nblow = l.second.first;
          int nbsegmentnumber = l.second.second;
          auto& nbsegment = segments[nbsegmentnumber];
          if (mylow==passheight && nblow<mylow && nbsegment.drainto!=mysegmentnumber)
            {
              // std::cout << " from segment " << it2->second << " to segment " << sorted_links[0].second.second
              //           << " my low=" << data[s.seedindex]
              //           << " passheight=" << s.links[sorted_links[0].second.second].passheight
              //           << " nb low=" << data[segments[sorted_links[0].second.second].seedindex]
              //           << " my accumulation=" << s.accumulation
              //           << " nb accumulation=" << segments[sorted_links[0].second.second].accumulation
              //           << std::endl;
              counter1++;
              mysegment.drainto = nbsegmentnumber;
              break;
            }
        }
      if (mysegment.drainto!=-1) continue; // we assigned a drainage

      // criterion 2:
      for (auto& l : sorted_links)
        {
          int passheight = l.first;
          int nblow = l.second.first;
          int nbsegmentnumber = l.second.second;
          auto& nbsegment = segments[nbsegmentnumber];
          if (mylow<passheight && nblow<mylow && nbsegment.drainto!=mysegmentnumber)
            {
              // std::cout << " from segment " << it2->second << " to segment " << sorted_links[0].second.second
              //           << " my low=" << data[s.seedindex]
              //           << " passheight=" << s.links[sorted_links[0].second.second].passheight
              //           << " nb low=" << data[segments[sorted_links[0].second.second].seedindex]
              //           << " my accumulation=" << s.accumulation
              //           << " nb accumulation=" << segments[sorted_links[0].second.second].accumulation
              //           << std::endl;
              counter2++;
              mysegment.drainto = nbsegmentnumber;
              break;
            }
        }
      if (mysegment.drainto!=-1) continue; // we assigned a drainage
    }
  std::cout << "low hanging criterion 1 drains " << counter1 << " segments" << std::endl; 
  std::cout << "low hanging criterion 2 drains " << counter2 << " segments" << std::endl;

  // update accumulation on pixel level
  std::cout << "update accumulation" << std::endl;
  for (auto it2=xs.rbegin(); it2!=xs.rend(); ++it2)
    {
      int mylow = it2->first;
      int mysegmentnumber = it2->second;
      auto& mysegment = segments[mysegmentnumber]; // this is our segment

      // does this segment drain to another segment?
      if (mysegment.drainto==-1)
        continue;

      // get links
      auto nbsegmentnumber = mysegment.drainto;
      auto& nbsegment = segments[mysegment.drainto];
      auto mylinkit = mysegment.links.find(nbsegmentnumber);
      if (mylinkit==mysegment.links.end())
        {
          std::cout << "Uuups cannot find that neighbor" << std::endl;
          continue;
        }
      auto nblinkit = nbsegment.links.find(mysegmentnumber);
      if (nblinkit==nbsegment.links.end())
        {
          std::cout << "Uuups cannot find that neighbor" << std::endl;
          continue;
        }

      // std::cout << " from " << mysegmentnumber << " to " << nbsegmentnumber
      //           << " fromheight=" << data[mysegment.seedindex]
      //           << " mypass=" << data[mylinkit->second.index]
      //           << " nbpass=" << data[nblinkit->second.index]
      //           << " toheight=" << data[nbsegment.seedindex]
      //           << " fromtype=" << mysegment.type << " totype=" << nbsegment.type
      //           << std::endl;
      // update accumulation on pixel level
      // mysegment(.seedindex) -> passheight
      if (mysegment.type==basin)
        {
          // modify path from pit to this side of the passheight
          auto newheight = data[mysegment.seedindex];
          std::vector<int> path_to_pass;
          path_to_pass.push_back(mylinkit->second.index);
          // std::cout << "X=" << mylinkit->second.index << ": " << mysegment.seedindex << " ";
          while (path_to_pass.back()!=mysegment.seedindex)
            {
              // std::cout << drainto[path_to_pass.back()] << " ";
              path_to_pass.push_back(drainto[path_to_pass.back()]);
            }
          // std::cout << std::endl;
          for (int k=path_to_pass.size()-2; k>=0; k--)
            {
              int fromindex = path_to_pass[k+1];
              int toindex = path_to_pass[k];
              accumulation[fromindex] -= accumulation[toindex];
              accumulation[toindex] += accumulation[fromindex];
              drainto[fromindex] = toindex;
              drainto[toindex] = -1;
              data[toindex] = newheight;
            }
          if (nbsegment.type==plain)
            nbsegment.accumulation += accumulation[mylinkit->second.index]; // accumulate here, because we do not know yet  where nb drains
        }
      else if (mysegment.type==plain)
        {
          // this is the pixel where the plain drains;
          // this pixel should have accumulation 1 at this time; now it gets the whole segment
          accumulation[mylinkit->second.index] = mysegment.accumulation + mysegment.size; // now we know where it drains
          if (nbsegment.type==plain)
            nbsegment.accumulation += accumulation[mylinkit->second.index]; // accumulate here, because we do not know yet where nb drains
        }
      // passheight -> nbsegment(.seedindex)
      if (nbsegment.type==basin)
        {
          // modify path from neighbors side of the pass height to the neighbors pit
          auto mypassheightindex = mylinkit->second.index;
          auto nbpassheightindex = nblinkit->second.index;
          auto newheight = data[mypassheightindex]; // this is the height of my pit to the pass height
          drainto[mypassheightindex] = nbpassheightindex;
          accumulation[nbpassheightindex] += accumulation[mypassheightindex];
          data[nbpassheightindex] = std::min(data[nbpassheightindex],newheight);
          std::vector<int> path_to_pit;
          path_to_pit.push_back(nbpassheightindex); // the start point
          // std::cout << "Z=" << nbsegment.seedindex << ": " << nbpassheightindex << " ";
          while (path_to_pit.back()!=nbsegment.seedindex)
            {
              // std::cout << drainto[path_to_pit.back()] << " ";
              path_to_pit.push_back(drainto[path_to_pit.back()]);
            }
          // std::cout << std::endl;
          for (int k=1; k<path_to_pit.size(); k++)
            {
              int fromindex = path_to_pit[k-1];
              int toindex = path_to_pit[k];
              accumulation[toindex] += accumulation[mypassheightindex];
              data[toindex] = std::min(data[toindex],newheight);
            }
        }
      else if (nbsegment.type==plain)
        {
          // we do not know yet where this segment drains
          // the mass accumulating to this segment is in the segments accumulation
        }
    }
  
  // show segment average height
  std::cout << "result" << std::endl;
  for (int i=0; i<N; i++)
    result[i] = accumulation[i];
  return result_image;

  // pits, plains and passes
  std::vector<int> ppp(N,0); // interior
  std::cout << "pits, plains and passes" << std::endl;
  for (int j=0; j<n; j++)
    for (int i=0; i<m; i++)
      {
        int index = j*m+i;
        if (drainto[index]<0 && accumulation[index]>1)
          {
            ppp[index] = -1; // a pit
            continue;
          }
        if (drainto[index]<0 && accumulation[index]==1)
          {
            ppp[index] = 1; // plain
            continue;
          }
        // now we are at a pixel which drains somewhere
        int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0));
        for (int k=0; k<neighbors[c].size(); k++)
          {
            int index_neighbor=index+neighbors[c][k];
            if (segment[index_neighbor]!=segment[index])
              ppp[index] = 2; //data[index]-data[segments[segment[index]].seedindex];
          }
      }
  

  // sort planar segments by size in order to find lakes
  std::vector<X> findlakes;
  for (int i=0; i<segments.size(); i++)
    if (segments[i].type==plain)
      findlakes.push_back(X(segments[i].size,i));
  std::sort(findlakes.begin(),findlakes.end());
  for (auto it=findlakes.rbegin(); it!=findlakes.rend(); ++it)
    std::cout << it->first << " " << it->second << std::endl;
  
  std::cout << "finished" << std::endl;
  result = ppp;
  return result_image;
}


template<typename T>
RasterDataSet<T> watershed (const RasterDataSet<T>& image)
{
  // ensure pixel type is int
  if (std::is_same<T,int>::value==false)
    {
      std::cout << "pixel type must be int" << std::endl;
      exit(1);
    }
    
  // prepare input
  const std::vector<T>& data=image.data(); // access to raw image data
  int N=data.size(); // number of pixels in image
  int m=image.sizeLong(); // pixels per row
  int n=image.sizeLat(); // pixels per column
  std::cout << "N=" << N << " m=" << m << " n=" << n << std::endl;

  /* prepare cases
   *  7 6 8
   *  1 0 2
   *  4 3 5
   */
  std::vector<std::vector<int>> neighbors ={{-m-1,-m,-m+1,-1,1,m-1,m,m+1}, // 0
                                            {-m,-m+1,1,m,m+1}, // 1
                                            {-m-1,-m,-1,m-1,m}, // 2
                                            {-1,1,m-1,m,m+1}, // 3
                                            {1,m,m+1}, // 4
                                            {-1,m-1,m}, // 5
                                            {-m-1,-m,-m+1,-1,1}, // 6
                                            {-m,-m+1,1}, // 7
                                            {-m-1,-m,-1}}; // 8
  std::vector<std::vector<int>> neighbors5 ={{-m,-1,1,m}, // 0
                                             {-m,1,m}, // 1
                                             {-m,-1,m}, // 2
                                             {-1,1,m}, // 3
                                             {1,m}, // 4
                                             {-1,m}, // 5
                                             {-m,-1,1}, // 6
                                             {-m,1}, // 7
                                             {-m,-1}}; // 8

  // prepare output
  RasterDataSet<T> result_image(image);
  std::vector<T>& result = result_image.data();
  T fillvalue = -9999;
  for (auto& pixel : result) pixel=fillvalue;

  // find all local minima
  // a local minimum is a pixel which has no neighbor with smaller height
  std::cout << "find minima" << std::endl;
  using X = std::pair<int,int>;
  std::priority_queue<X, std::vector<X>, std::greater<X>> sinks;
  for (int j=0; j<n; j++)
    for (int i=0; i<m; i++)
      {
        int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0));
        int index = j*m+i;
        T value = data[index];
        int down=0;
        for (int k=0; k<neighbors[c].size(); k++)
          if (data[index+neighbors[c][k]]<value) down++;
        if (down==0)
          {
            //std::cout << "push: value=" << value << " index=" << index << std::endl;
            sinks.push(X(value,index));
          }
      }
  std::cout << "minima=" << sinks.size() << "=" << (int)round(100.0*sinks.size()/N) << "%"
            << std::endl;

  
  // now classify pixels and build segments
  std::cout << "find initial segments" << std::endl;
  T segments=0;
  int classified=0;
  std::vector<int> segment_lowestindex; // the seed index for the segment
  std::vector<int> segment_highestindex; // index of pixel with highest value
  while (!sinks.empty())
    {
      X x = sinks.top();
      int value = x.first;
      int seedindex = x.second;
      // std::cout << "pop: value=" << value << " index=" << seedindex << " result=" << result[seedindex] << std::endl;
      sinks.pop();
      if (result[seedindex]!=fillvalue) continue;

      // pixel is not yet classified; start a new segment
      segment_lowestindex.push_back(seedindex);
      std::queue<int> q;
      q.push(seedindex);
      T maxheight=-100000;
      int maxindex;
      while (!q.empty())
        {
          // get next pixel from queue
          int index = q.front(); q.pop();
          // if pixel is already classified then skip it
          if (result[index]!=fillvalue) continue;
          // otherwise assign it to the current segment
          result[index] = segments;
          classified++;
          if (data[index]>maxheight)
            {
              maxheight=data[index];
              maxindex=index;
            }
          // now go to neighbors which are not lower
          int i=index%m;
          int j=index/m;
          int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0));
          for (int k=0; k<neighbors[c].size(); k++)
            {
              int index_neighbor=index+neighbors[c][k];
              if (data[index_neighbor]>=data[index] && result[index_neighbor]==fillvalue)
                q.push(index_neighbor);
            }
        }
      // a segment has been finished
      segment_highestindex.push_back(maxindex);
      segments++;
    }
  std::cout << "segments=" << segments << " classified=" << classified << " unclassified=" << N-classified << std::endl;

  // build histogram
  if (false)
    {
      std::cout << "compute histogram" << std::endl;
      int delta=10;
      std::vector<int> count(5000/delta+1,0);
      for (int i=0; i<segments; i++)
        {
          int difference=data[segment_highestindex[i]]-data[segment_lowestindex[i]];
          int bin=difference/delta;
          count[bin]++;
        }
      for (int i=0; i<count.size(); i++)
        std::cout << "bin " << i << " from " << i*delta << " to " << (i+1)*delta << " counts " << count[i] << std::endl;
    }
  
  // merge segments
  // build segment graph
  std::cout << "construct segment graph" << std::endl;
  std::vector<std::set<int>> segment_neighbors(segments); // list of all neighboring segments
  std::vector<std::map<int,int>> pass_height(segments); // height of lowest pixel between two segments
  std::vector<int> segment_size(segments); // number of pixels in segment
  std::vector<bool> visited(N,false);
  for (auto seedindex : segment_lowestindex)
    {
      // visit all pixels of that segment and determine neighbors
      std::queue<int> q;
      q.push(seedindex);
      int mysegment = result[seedindex];
      int count = 0;
      // std::cout << "start segment " << mysegment << " with seedindex " << seedindex << std::endl;
      while (!q.empty())
        {
          // get next pixel from queue
          int index = q.front(); q.pop();
          if (visited[index]) continue;
          visited[index] = true;
          count++;
          // now visit all neighbors
          int i=index%m;
          int j=index/m;
          int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0));
          for (int k=0; k<neighbors[c].size(); k++)
            {
              int index_neighbor = index+neighbors[c][k];
              int neighborsegment = result[index_neighbor];
              // std::cout << "index=" << index << " neighbor=" << index_neighbor << " mysegment=" << mysegment
              //           << " hissegment=" << neighborsegment << " visited=" << visited[index_neighbor] << std::endl;
              if (neighborsegment==mysegment)
                {
                  // neighbor is in the same segment
                  if (!visited[index_neighbor])
                    q.push(index_neighbor);
                }
              else
                {
                  // neighbor is in different segment; so we have a neighbor
                  segment_neighbors[mysegment].insert(neighborsegment);
                  auto it = pass_height[mysegment].find(neighborsegment);
                  if (it==pass_height[mysegment].end())
                    pass_height[mysegment][neighborsegment] = std::max(data[index],data[index_neighbor]);
                  else
                    it->second = std::max(it->second,std::max(data[index],data[index_neighbor]));
                }
            }
        }
      segment_size[mysegment] = count;
      // std::cout << "finished segment " << mysegment << " with " << count << " pixels" << std::endl;
    }
  for (int i=0; i<visited.size(); i++) visited[i]=false;
  // consistency check
  for (int thissegment=0; thissegment<segments; thissegment++)
    if (segment_size[thissegment]>0)
      {
        for (auto neighborsegment : segment_neighbors[thissegment])
          {
            if (segment_neighbors[neighborsegment].count(thissegment)!=1)
              std::cout << neighborsegment << " is neighbor of " << thissegment << " but not vice versa" << std::endl;
            if (pass_height[thissegment].find(neighborsegment)==pass_height[thissegment].end())
              std::cout << thissegment << " segment_neighbors and pass_height inconsistent" << std::endl;
          }
        for (auto it=pass_height[thissegment].begin(); it!=pass_height[thissegment].end(); ++it)
          {
            auto itreverse = pass_height[it->first].find(thissegment);
            if (itreverse==pass_height[it->first].end())
              std::cout << it->first << " is neighbor of " << thissegment << " in pass_height but not vice versa" << std::endl;
            else
              if (it->second!=itreverse->second)
                std::cout << "pass height not symmetric " << it->first << " " << itreverse->first << " "
                          << it->second << " " << itreverse->second << std::endl;
          }
      }
  
  // merge segments based on pass_height thresholds
  std::vector<int> thresholds = {10,10,10,10,25,25,25,50,50,50,75,75,75,75,75,75,75};
  for (int threshold : thresholds)
    {
      int merged_segments=0;
      for (int thissegment=0; thissegment<segments; thissegment++)
        if (segment_size[thissegment]>0)
          {
            // try to merge this segment to a neighbor
            // this neighbor may be before or after us in the list
            int thissegment_low = data[segment_lowestindex[thissegment]];
            int merge_with_segment=-1;
            for (auto neighborsegment : segment_neighbors[thissegment])
              {
                int neighborsegment_low = data[segment_lowestindex[neighborsegment]];
                int passheight = pass_height[thissegment][neighborsegment];
                if (neighborsegment_low<=thissegment_low && passheight-thissegment_low<threshold)
                  {
                    // we can be merged with neighbor
                    merge_with_segment = neighborsegment;
                    break;
                  }
              }

            if (merge_with_segment<0)
              continue; // nothing to merge

            // we will merge this segment
            merged_segments++;
            
            // merge this segment with a neighboring segment
            // delete access to this segment from all neighbors lists and add the merge_with
            //std::cout << "merging segment " << thissegment << " to segment " << merge_with_segment << std::endl;
            for (auto neighborsegment : segment_neighbors[thissegment])
              {
                // it is not the neighbor I will merge with
                // but this neighbor will be a neighbor to merge_with then
                // modify neighbor set in neighbor
                segment_neighbors[neighborsegment].erase(thissegment); // erase me
                auto it = pass_height[neighborsegment].find(thissegment);
                if (it==pass_height[neighborsegment].end())
                  std::cout << " ups should not happen" << std::endl;
                auto old_height=it->second;
                pass_height[neighborsegment].erase(it); // erase me

                // the rest is only needed if neighbor is not the segment we merge with
                if (neighborsegment==merge_with_segment)
                  continue;
                      
                // make new neighbors
                segment_neighbors[neighborsegment].insert(merge_with_segment); // add new neighbor (possibly)
                segment_neighbors[merge_with_segment].insert(neighborsegment); // add new neighbor (possibly)

                // modify pass height in neighbors
                auto it2 = pass_height[neighborsegment].find(merge_with_segment); // find new neighbor
                if (it2==pass_height[neighborsegment].end())
                  pass_height[neighborsegment][merge_with_segment] = old_height; // add new neighbor if not found
                else
                  it2->second = std::min(old_height,it->second); // compute new ridge height
                it2 = pass_height[merge_with_segment].find(neighborsegment); // find the other direction
                if (it2==pass_height[merge_with_segment].end())
                  pass_height[merge_with_segment][neighborsegment] = old_height; // add new neighbor if not found
                else
                  it2->second = std::min(old_height,it->second); // compute new ridge height
              }
        
            // delete neighbors of this segment
            segment_neighbors[thissegment].clear(); // no neighbors anymore
            pass_height[thissegment].clear(); // same same

            // adjust segment sizes
            segment_size[merge_with_segment] += segment_size[thissegment];
            segment_size[thissegment] = 0;

            // now the only thing left is to change the result image
            // we can use the result as visited indicator
            int seedindex = segment_lowestindex[thissegment];
            std::queue<int> q;
            q.push(seedindex);
            while (!q.empty())
              {
                // get next pixel from queue
                int index = q.front(); q.pop();
                if (result[index]!=thissegment) continue;
                result[index] = merge_with_segment; // this is the new segment
                // now visit all neighbors
                int i=index%m;
                int j=index/m;
                int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0));
                for (int k=0; k<neighbors[c].size(); k++)
                  {
                    int index_neighbor = index+neighbors[c][k];
                    int neighborsegment = result[index_neighbor];
                    if (neighborsegment==thissegment)
                      q.push(index_neighbor);
                  }
              }
            //std::cout << "merged segment " << thissegment << " to segment " << merge_with_segment << std::endl;
          }
      std::cout << "merged " << merged_segments << " segments at threshold " << threshold << std::endl;
      // check consistency
      for (int thissegment=0; thissegment<segments; thissegment++)
        if (segment_size[thissegment]>0)
          {
            for (auto neighborsegment : segment_neighbors[thissegment])
              {
                if (segment_neighbors[neighborsegment].count(thissegment)!=1)
                  std::cout << neighborsegment << " is neighbor of " << thissegment << " but not vice versa" << std::endl;
                if (pass_height[thissegment].find(neighborsegment)==pass_height[thissegment].end())
                  std::cout << thissegment << " segment_neighbors and pass_height inconsistent" << std::endl;
              }
            for (auto it=pass_height[thissegment].begin(); it!=pass_height[thissegment].end(); ++it)
              {
                auto itreverse = pass_height[it->first].find(thissegment);
                if (itreverse==pass_height[it->first].end())
                  std::cout << it->first << " is neighbor of " << thissegment << " in pass_height but not vice versa" << std::endl;
                else
                  if (it->second!=itreverse->second)
                    std::cout << "pass height not symmetric " << it->first << " " << itreverse->first << " "
                              << it->second << " " << itreverse->second << std::endl;
              }
          }
    }

  // check pass height
  if (false)
    for (int thissegment=0; thissegment<segments; thissegment++)
      if (segment_size[thissegment]>0)
        for (auto it=pass_height[thissegment].begin(); it!=pass_height[thissegment].end(); ++it)
          if (it->second-data[segment_lowestindex[thissegment]]<75)
            std::cout << "pass from " << thissegment << " to " << it->first << " is "
                      << it->second-data[segment_lowestindex[thissegment]] << std::endl;

  
  // clean up segment numbers
  std::vector<int> new_segment_number(segments);
  int count = 0;
  for (int thissegment=0; thissegment<segments; thissegment++)
    if (segment_size[thissegment]>0)
      new_segment_number[thissegment] = count++;
  for (int i=0; i<N; i++) result[i] = new_segment_number[result[i]];
  
  std::cout << "finished with " << count << " segments" << std::endl;
  return result_image;
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

    // open ini file
    Dune::ParameterTree ptree;
    Dune::ParameterTreeParser ptreeparser;
    // ptreeparser.readINITree("readtiff.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);

    // read precipitation data from netcdf file
    if (false)
      {
        NetCDFFile netcdf("rr_0.1deg_day_2020_grid_ensmean.nc","/Users/peterbastian/Data/incoming/");
        std::cout << "global attributes" << std::endl;
        std::cout << netcdf.get_string_attribute(0) << std::endl;
        std::cout << netcdf.get_string_attribute(1) << std::endl;
        std::cout << netcdf.get_string_attribute(2) << std::endl;
        std::cout << netcdf.get_string_attribute(3) << std::endl;
        std::cout << "units" << std::endl;
        std::cout << netcdf.get_string_attribute(0,2) << std::endl;
        std::cout << netcdf.get_string_attribute(1,2) << std::endl;
        std::cout << netcdf.get_string_attribute(2,2) << std::endl;
        std::cout << "variable 3 attributes" << std::endl;
        std::cout << netcdf.get_string_attribute(3,0) << std::endl;
        std::cout << netcdf.get_string_attribute(3,1) << std::endl;
        std::cout << netcdf.get_string_attribute(3,2) << std::endl;
        std::cout << netcdf.get_attribute<float>(3,3) << std::endl;
        std::cout << netcdf.get_attribute<float>(3,4) << std::endl;
        std::cout << netcdf.get_attribute<short>(3,5) << std::endl;
        std::cout << netcdf.get_attribute<short>(3,6) << std::endl;

        auto v_long = netcdf.get_variable<double>(0);
        print_vector_overview(v_long);
        auto v_lat = netcdf.get_variable<double>(1);
        print_vector_overview(v_lat);
        auto v_time = netcdf.get_variable<double>(2);
        print_vector_overview(v_time);
        auto FillValue = netcdf.get_attribute<short>(3,5);
        std::cout << "FillValue is " << FillValue << std::endl;

        const int dim=2;
        Dune::FieldVector<double,dim> L;
        L[0] = netcdf.dimension_length(0);
        L[1] = netcdf.dimension_length(1);
        std::array<int,dim> N;
        N[0] = netcdf.dimension_length(0);;
        N[1] = netcdf.dimension_length(1);;
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        auto gridp = std::make_shared<Grid>(L,N,std::bitset<dim>(0ULL),1);
        typedef Grid::LeafGridView GV;
        GV gv=gridp->leafGridView();

        // write a grid file
        std::string filename("rr");
        struct stat st;
        if( stat( filename.c_str(), &st ) != 0 )
          {
            int stat = 0;
            stat = mkdir(filename.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
            if( stat != 0 && stat != -1)
              std::cout << "Error: Cannot create directory "
                        << filename << std::endl;
          }
        typedef Dune::VTKWriter<GV> VTKWRITER;
        VTKWRITER vtkwriter(gv,Dune::VTK::conforming);
        typedef Dune::VTKSequenceWriter<GV> VTKSEQUENCEWRITER;
        VTKSEQUENCEWRITER vtkSequenceWriter(std::make_shared<VTKWRITER>(vtkwriter),filename,filename,"");
        auto rr = netcdf.get_variable<short>(3,0);
        vtkSequenceWriter.addCellData(rr,"rainfall [mm]");
        float offset = netcdf.get_attribute<float>(3,3);
        float scale_factor = netcdf.get_attribute<float>(3,4);
        for (int k=0; k<netcdf.dimension_length(2); k++)
          {
            auto rrk = netcdf.get_variable<short>(3,k);
            RasterDataSet<short>(v_long[0],v_lat[0],v_long[1]-v_long[0],v_lat[1]-v_lat[0],netcdf.dimension_length(0),netcdf.dimension_length(1),rrk,1);
            for (int i=0; i<rr.size(); i++)
              rr[i] = offset + scale_factor*std::max((short)0,rrk[i]);
            vtkSequenceWriter.write(k,Dune::VTK::appendedraw);
          }
      }

    if (argc==2)
      {
        std::string file(argv[1]);
        GeoTIFFInfo image(file,"",1);
        // const int dim=2;
        // Dune::FieldVector<double,dim> L;
        // L[0] = 1.0;
        // L[1] = 1.0;
        // std::array<int,dim> N;
        // N[0] = image.sizeLong();
        // N[1] = image.sizeLat();
        // typedef Dune::YaspGrid<dim> Grid;
        // typedef Grid::ctype DF;
        // auto gridp = std::make_shared<Grid>(L,N,std::bitset<dim>(0ULL),1);
        return 0;
      }

    // Create BWat90m data set
    if (true)
      {
        // read and access tiff files for Baden-Württemberg
        // GeoTIFFImage<GInt16> srtm_38_03("srtm_38_03.tif","/Users/peterbastian/Data/BWat90m/",1,2);
        // GeoTIFFImage<GInt16> srtm_39_03("srtm_39_03.tif","/Users/peterbastian/Data/BWat90m/",1,2);
        GeoTIFFImage<GByte> land_cover("bwlandcover.tif","/Users/peterbastian/Data/BWat90m/",1,2);

        GeoTIFFImage<GInt16> hs_con_n45e005("hs_con_n45e005.tif","/Users/peterbastian/Data/HydroSHEDS/eu_con_3s_zip_grid/n45e005_con_grid/n45e005_con/n45e005_con/",1,2);
        GeoTIFFImage<GInt16> hs_con_n45e010("hs_con_n45e010.tif","/Users/peterbastian/Data/HydroSHEDS/eu_con_3s_zip_grid/n45e010_con_grid/n45e010_con/n45e010_con/",1,2);

        // make YaspGrid for a 4 by 3 degree simulation area with resolution 1/1200 degrees
        double dx=std::abs(hs_con_n45e005.dLong());
        double dy=std::abs(hs_con_n45e005.dLat());
        const int dim=2;
        Dune::FieldVector<double,dim> L;
        L[0] = 4.0;
        L[1] = 4.0;
        std::array<int,dim> N;
        N[0] = 4*1200;
        N[1] = 4*1200;
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        auto gridp = std::make_shared<Grid>(L,N,std::bitset<dim>(0ULL),1);
        typedef Grid::LeafGridView GV;
        GV gv=gridp->leafGridView();


        // make a P0 grid function
        typedef typename GV::Grid::ctype DF; // type for ccordinates
        typedef double RF;                   // type for computations
        typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
        FEM fem(Dune::GeometryTypes::cube(dim));
        typedef Dune::PDELab::NoConstraints CON;
        typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
        typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
        GFS gfs(gv,fem);
        gfs.name("Q0");

        using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
        using ZDGF = Dune::PDELab::DiscreteGridFunction<GFS,Z>;
        using VTKF = Dune::PDELab::VTKGridFunctionAdapter<ZDGF>;

        Z r1(gfs);
        ZDGF r1dgf(gfs,r1);
        Z r2(gfs);
        ZDGF r2dgf(gfs,r2);
        Z r3(gfs);
        ZDGF r3dgf(gfs,r3);
        Z r4(gfs);
        ZDGF r4dgf(gfs,r4);
        Z r5(gfs);
        ZDGF r5dgf(gfs,r5);

        // now make raster canvas in cell-centered mode
        RasterDataSet<int> elevation(7.0+0.5*dx,46.0+0.5*dy,dx,dy,N[0],N[1],0.0,1);
        elevation.paste(hs_con_n45e005);
        elevation.paste(hs_con_n45e010);
        RasterDataSet<char> landcover(7.0+0.5*dx,46.0+0.5*dy,dx,dy,N[0],N[1],-1.0,1);
        landcover.paste(land_cover);

        // auto glambda = [](const auto& x){auto x0=M_PI*x[0]/4.0; auto x1=M_PI*x[1]/4.0;
        //                                  return sin(9*x0)*cos(10*x1);};
        // auto g = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda);
        // Dune::PDELab::interpolate(g,gfs,z);
        // for (size_t i=0; i<elevation.data().size(); i++)
        //   elevation.data()[i] = floor(1001+1000*Dune::PDELab::Backend::native(z)[i][0]);

        // auto result=watershed(elevation);
        elevation = outlier(elevation,2,2000);
        auto result2=geometric_average(elevation,3);
        auto result3=variance(elevation,3);
        auto result4=median(elevation,3);

        for (size_t i=0; i<elevation.data().size(); i++)
          Dune::PDELab::Backend::native(r1)[i][0] = elevation.data()[i];
        for (size_t i=0; i<elevation.data().size(); i++)
          Dune::PDELab::Backend::native(r2)[i][0] = result2.data()[i];
        for (size_t i=0; i<elevation.data().size(); i++)
          Dune::PDELab::Backend::native(r3)[i][0] = result3.data()[i];
        for (size_t i=0; i<elevation.data().size(); i++)
          Dune::PDELab::Backend::native(r4)[i][0] = result4.data()[i];
        for (size_t i=0; i<elevation.data().size(); i++)
          Dune::PDELab::Backend::native(r5)[i][0] = landcover.data()[i];

        using Writer = Dune::VtkImageDataWriter<GV>;
        Writer vtkWriter(gv, Dune::Vtk::FormatTypes::COMPRESSED, Dune::Vtk::DataTypes::FLOAT32);
        vtkWriter.addCellData(std::make_shared<VTKF>(r1dgf,"elevation"));
        vtkWriter.addCellData(std::make_shared<VTKF>(r2dgf,"average"));
        vtkWriter.addCellData(std::make_shared<VTKF>(r3dgf,"variance"));
        vtkWriter.addCellData(std::make_shared<VTKF>(r4dgf,"median"));
        vtkWriter.addCellData(std::make_shared<VTKF>(r5dgf,"landcover"));
        vtkWriter.write("BWat90m.vti");

        return 0;
      }

    // Heidelbergat25m data set
    if (false)
      {
        // GeoTIFFImage<float> image_03_02("EUD_CP-DEMS_4500025000-AA_WGS84_03_02.tif","/Volumes/Samsung_T5/Data/tiles/",1,2);
        // GeoTIFFImage<float> image_03_03("EUD_CP-DEMS_4500025000-AA_WGS84_03_03.tif","/Volumes/Samsung_T5/Data/tiles/",1,2);
        // GeoTIFFImage<float> image_03_04("EUD_CP-DEMS_4500025000-AA_WGS84_03_04.tif","/Volumes/Samsung_T5/Data/tiles/",1,2);
        GeoTIFFImage<float> image_02_02("EUD_CP-DEMS_4500025000-AA_WGS84_02_02.tif","/Volumes/Samsung_T5/Data/tiles/",1,2);
        GeoTIFFImage<float> image_02_03("EUD_CP-DEMS_4500025000-AA_WGS84_02_03.tif","/Volumes/Samsung_T5/Data/tiles/",1,2);
        GeoTIFFImage<float> image_02_04("EUD_CP-DEMS_4500025000-AA_WGS84_02_04.tif","/Volumes/Samsung_T5/Data/tiles/",1,2);
        GeoTIFFImage<float> image_01_02("EUD_CP-DEMS_4500025000-AA_WGS84_01_02.tif","/Volumes/Samsung_T5/Data/tiles/",1,2);
        GeoTIFFImage<float> image_01_03("EUD_CP-DEMS_4500025000-AA_WGS84_01_03.tif","/Volumes/Samsung_T5/Data/tiles/",1,2);
        GeoTIFFImage<float> image_01_04("EUD_CP-DEMS_4500025000-AA_WGS84_01_04.tif","/Volumes/Samsung_T5/Data/tiles/",1,2);
        GeoTIFFImage<float> image_06_03("EUD_CP-DEMS_4500035000-AA_WGS84_06_03_resampled.tif","/Volumes/Samsung_T5/Data/tiles/",1,2);
        GeoTIFFImage<float> image_06_04("EUD_CP-DEMS_4500035000-AA_WGS84_06_04_resampled.tif","/Volumes/Samsung_T5/Data/tiles/",1,2);        

        // make YaspGrid for a 4 by 3 degree simulation area with resolution 1/1200 degrees
        // Heidelberg = (8.66,49.4), compute position of Heidelberg on the raster
        double dx=std::abs(image_02_02.dLong());
        double dy=std::abs(image_02_02.dLat());
        double ox = image_02_02.originLong();
        double oy = image_02_02.originLat();
        double hdx = ox+round((8.66-ox)/dx)*dx;
        double hdy = oy+round((49.4-oy)/dy)*dy;
        const int dim=2;
        Dune::FieldVector<double,dim> L;
        L[0] = 3000*dx;
        L[1] = 2000*dy;
        std::array<int,dim> N;
        N[0] = 3000;
        N[1] = 2000;
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        auto gridp = std::make_shared<Grid>(L,N,std::bitset<dim>(0ULL),1);

        // now make raster canvas in cell-centered mode
        RasterDataSet<float> canvas_elev(hdx-1000*dx,hdy-1000*dy,dx,dy,N[0],N[1],0.0,1);
        canvas_elev.paste(image_01_02);
        canvas_elev.paste(image_01_03);
        canvas_elev.paste(image_01_04);
        canvas_elev.paste(image_02_02);
        canvas_elev.paste(image_02_03);
        canvas_elev.paste(image_02_04);
        canvas_elev.paste(image_06_03);
        canvas_elev.paste(image_06_04);

        // write a grid file
        typedef Grid::LeafGridView GV;
        GV gv=gridp->leafGridView();
        typedef Dune::VTKWriter<GV> VTKWRITER;
        VTKWRITER vtkwriter(gv,Dune::VTK::conforming);
        vtkwriter.addCellData(canvas_elev.data(),"elevation");
        vtkwriter.write("HDat25m",Dune::VTK::appendedraw);
        
        return 0;
      }
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (GeoTIFFReaderException &e){
    std::cerr << "geotiff error: " << e.what() << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
