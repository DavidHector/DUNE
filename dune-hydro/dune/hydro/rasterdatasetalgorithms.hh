#ifndef Udune_hydro_rasterdataset_algorithms_HH
#define Udune_hydro_rasterdataset_algorithms_HH

#include <exception>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include"rasterdataset.hh"


// a function for outlier correction in bathymmetry
template<typename T>
RasterDataSet<T> outlier (const RasterDataSet<T>& image_in, int diameter, int threshold)
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
            std::cout << "outlier: i=" << i << " j=" << j << " value=" << in[j*m+i] << " median=" << median << std::endl;
            out[j*m+i] = median;
          }
        else
          out[j*m+i] = in[j*m+i];
      }
  return image_out;
}

// computes number of pixels draining to a pixel
template<typename T>
RasterDataSet<int> drainage_accumulation (const RasterDataSet<T>& image)
{
  // analyse size etc of input image
  auto& data=image.data(); // access to raw image data
  int N=data.size(); // number of pixels in image
  int m=image.sizeLong(); // pixels per row
  int n=image.sizeLat(); // pixels per column
  std::cout << "accumulate: N=" << N << " m=" << m << " n=" << n << std::endl;

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

  // allocate output image
  RasterDataSet<int> result_image(image.originLong(),image.originLat(),image.dLong(),image.dLat(),
                                  image.sizeLong(),image.sizeLat(),0.0,1);
  auto& result = result_image.data();
  
  // compute accumulation and drainto fields
  std::vector<int> drainto(N,-1); // a pixel drains to the lowest neighbor that is strictly lower than itself;
                                  // -1 means there is none
  std::vector<int> accumulation(N,1); // number of other pixels draining to the given pixel + 1 (for itself)

  // sort pixels by height
  std::cout << "sort by height" << std::endl;
  using X = std::pair<int,int>;
  std::vector<X> x; // pixels sorted by height 
  for (int j=0; j<n; j++) // loop over rows
    for (int i=0; i<m; i++) // loop over columns
      {
        int index = j*m+i;
        T value = data[index];
        x.push_back(X(value,index));
      }
  std::sort(x.begin(),x.end());

  // now accumulate by proceeding from high to low!
  std::cout << "accumulate" << std::endl;
  for (auto it=x.rbegin(); it!=x.rend(); ++it) // this is important that we proceed from high to low
    {
      int height = it->first;
      int index = it->second;
      int i=index%m; // column
      int j=index/m; // row
      int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0)); // pixel boundary type
      int minheight=height;
      int mink=-1;
      // determine lowest lying neighbor
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
  // first statistics : count pits and plains only if they are *not* at the boundary
  int count_pits = 0;
  int count_plains = 0;
  for (int index=0; index<N; index++)
    {
      int i=index%m; // column
      int j=index/m; // row
      int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0)); // pixel boundary type
      if (c==0 && drainto[index]<0 && accumulation[index]>1) count_pits++;
      if (c==0 && drainto[index]<0 && accumulation[index]==1) count_plains++;
    }
  std::cout << count_pits << " pitpixels out of " << N << " pixels (" << (int)round(100.0*count_pits/N) << "%)" << std::endl;
  std::cout << count_plains << " plainpixels out of " << N << " pixels (" << (int)round(100.0*count_plains/N) << "%)" << std::endl;

  result = accumulation; // copy
  return result_image;
}

/* Given a direction map in ESRI format
 *
 * 32 64 128
 * 16  0   1
 *  8  4   2
 *
 * compute an image with distances from the outflow
 */
template<typename T>
RasterDataSet<int> distance_from_direction (const RasterDataSet<T>& image)
{
  // analyse size etc of input image
  auto& data=image.data(); // access to raw image data
  int N=data.size(); // number of pixels in image
  int m=image.sizeLong(); // pixels per row
  int n=image.sizeLat(); // pixels per column
  std::cout << "distance_from_direction: N=" << N << " m=" << m << " n=" << n << std::endl;

  /* prepare cases
   *  7 6 8
   *  1 0 2
   *  4 3 5
   */
  std::vector<std::vector<int>> neighbors ={
    {1,-m+1,-m,-m-1,-1,m-1,m,m+1}, // 0
    {1,-m+1,-m,m,m+1}, // 1
    {-m,-m-1,-1,m-1,m}, // 2
    {1,-1,m-1,m,m+1}, // 3
    {1,m,m+1}, // 4
    {-1,m-1,m}, // 5
    {1,-m+1,-m,-m-1,-1}, // 6
    {1,-m+1,-m}, // 7
    {-m,-m-1,-1}, // 8
  };
  std::vector<std::vector<int>> dir = {
    {1,2,4,8,16,32,64,128}, // 0
    {1,2,4,64,128}, // 1
    {4,8,16,32,64}, // 2
    {1,16,32,64,128}, // 3
    {1,64,128}, // 4
    {16,32,64}, // 5
    {1,2,4,8,16}, // 6
    {1,2,4}, // 7
    {4,8,16}, // 8
  };
  std::vector<int> backdir(129,0);
  backdir[1] = 16;
  backdir[2] = 32;
  backdir[4] = 64;
  backdir[8] = 128;
  backdir[16] = 1;
  backdir[32] = 2;
  backdir[64] = 4;
  backdir[128] = 8;
  std::vector<int> offset(129,0);
  offset[1] = 1;
  offset[2] = -m+1;
  offset[4] = -m;
  offset[8] = -m-1;
  offset[16] = -1;
  offset[32] = m-1;
  offset[64] = m;
  offset[128] = m+1;

  // allocate output image initialized with 1
  RasterDataSet<int> result_image(image.originLong(),image.originLat(),image.dLong(),image.dLat(),
                                  image.sizeLong(),image.sizeLat(),0,0);
  auto& distance = result_image.data();

  // we need to climb up the tree to compute the distance from the outlet
  std::vector<int> front; // pixels to treat in next round
  // determine roots, i.e. pixels not draining to another pixel
  for (int index=0; index<N; index++)
    {
      if (data[index]<=0) {front.push_back(index); continue;} // this is a sink
      int i=index%m; // column
      int j=index/m; // row
      if (i==0 && (data[index]>=8 && data[index]<=32))
        {
          // drains to outside pixel left
          front.push_back(index);
          continue;
        }
      if (i==m-1 && (data[index]<=2 || data[index]>=128))
        {
          // drains to outside pixel right
          front.push_back(index);
          continue;
        }
      if (j==0 && (data[index]>=2 && data[index]<=8))
        {
          // drains to outside pixel below
          front.push_back(index);
          continue;
        }
      if (j==n-1 && (data[index]>=32 && data[index]<=128))
        {
          // drains to outside pixel above
          front.push_back(index);
          continue;
        }
    }

  // climb up the tree and determine distance
  while (front.size()>0)
    {
      std::vector<int> newfront; // pixels to treat in next round
      for (auto index : front)
        {
          // get the current pixel
          int i=index%m; // column
          int j=index/m; // row
          int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0)); // pixel boundary type
          
          // loop through all neighbors of current pixel
          for (int k=0; k<neighbors[c].size(); k++)
            {
              int index_neighbor=index+neighbors[c][k];
              if (data[index_neighbor]==backdir[dir[c][k]])
                {
                  // this neighbor drains to the current pixel
                  if (distance[index_neighbor]>0)
                    {
                      std::cout << "ERROR pixel was already visited distance=" << distance[index_neighbor] << std::endl;
                                                                                                              exit(1);
                    }
                  distance[index_neighbor] = distance[index]+1;
                  newfront.push_back(index_neighbor);
                }
            }
        }
      front = newfront;
    }
  
  return result_image;
}

/* Given a direction map in ESRI format
 *
 * 32 64 128
 * 16  0   1
 *  8  4   2
 *
 * and a set of pixels compute an image of all pixels draining to the given set
 */
template<typename T>
RasterDataSet<int> upstream_from_direction (const RasterDataSet<T>& direction_image, std::vector<std::pair<int,int>> outlet)
{
  // analyse size etc of input image
  auto& direction_data=direction_image.data(); // access to raw image data
  int N=direction_data.size(); // number of pixels in image
  int m=direction_image.sizeLong(); // pixels per row
  int n=direction_image.sizeLat(); // pixels per column
  std::cout << "upstream_from_direction: N=" << N << " m=" << m << " n=" << n << std::endl;

  std::vector<std::vector<int>> neighbors ={
    {1,-m+1,-m,-m-1,-1,m-1,m,m+1}, // 0
    {1,-m+1,-m,m,m+1}, // 1
    {-m,-m-1,-1,m-1,m}, // 2
    {1,-1,m-1,m,m+1}, // 3
    {1,m,m+1}, // 4
    {-1,m-1,m}, // 5
    {1,-m+1,-m,-m-1,-1}, // 6
    {1,-m+1,-m}, // 7
    {-m,-m-1,-1}, // 8
  };
  std::vector<std::vector<int>> dir = {
    {1,2,4,8,16,32,64,128}, // 0
    {1,2,4,64,128}, // 1
    {4,8,16,32,64}, // 2
    {1,16,32,64,128}, // 3
    {1,64,128}, // 4
    {16,32,64}, // 5
    {1,2,4,8,16}, // 6
    {1,2,4}, // 7
    {4,8,16}, // 8
  };
  std::vector<int> backdir(129,0);
  backdir[1] = 16;
  backdir[2] = 32;
  backdir[4] = 64;
  backdir[8] = 128;
  backdir[16] = 1;
  backdir[32] = 2;
  backdir[64] = 4;
  backdir[128] = 8;

  // allocate output image initialized with zero
  RasterDataSet<int> result_image(direction_image.originLong(),direction_image.originLat(),
                                  direction_image.dLong(),direction_image.dLat(),
                                  direction_image.sizeLong(),direction_image.sizeLat(),0,0);
  auto& result_data = result_image.data();

  // seed front
  std::vector<int> front; // pixels to treat in next round
  for (auto pixel : outlet)
    {
      int i = pixel.first;
      int j = pixel.second;
      int index = j*m+i;
      if (result_data[index]==0)
        {
          result_data[index] = 1;
          front.push_back(index);
        }
    }
  std::cout << front.size() << " pixels in front" << std::endl;
  
  // breadth first graph search
  while (front.size()>0)
    {
      std::cout << front.size() << " pixels in front" << std::endl;
      std::vector<int> newfront; // pixels to treat in next round
      for (auto index : front)
        {
          // get the current pixel
          int i=index%m; // column
          int j=index/m; // row
          int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0)); // pixel boundary type
          
          // loop through all neighbors of current pixel
          for (int k=0; k<neighbors[c].size(); k++)
            {
              int index_neighbor=index+neighbors[c][k];
              if (result_data[index_neighbor]==0 && direction_data[index_neighbor]==backdir[dir[c][k]])
                {
                  result_data[index_neighbor] = 1;
                  newfront.push_back(index_neighbor);
                }
            }
        }
      front = newfront;
    }

  // return result
  return result_image;
}

// computes number of pixels draining to a pixel from a direction map in ESRI format
template<typename T>
RasterDataSet<int> accumulation_from_direction (const RasterDataSet<T>& image)
{
  // analyse size etc of input image
  auto& data=image.data(); // access to raw image data
  int N=data.size(); // number of pixels in image
  int m=image.sizeLong(); // pixels per row
  int n=image.sizeLat(); // pixels per column
  std::cout << "accumulation_from_direction: N=" << N << " m=" << m << " n=" << n << std::endl;

  /* prepare cases
   *  7 6 8
   *  1 0 2
   *  4 3 5
   */
  std::vector<std::vector<int>> neighbors ={
    {1,-m+1,-m,-m-1,-1,m-1,m,m+1}, // 0
    {1,-m+1,-m,m,m+1}, // 1
    {-m,-m-1,-1,m-1,m}, // 2
    {1,-1,m-1,m,m+1}, // 3
    {1,m,m+1}, // 4
    {-1,m-1,m}, // 5
    {1,-m+1,-m,-m-1,-1}, // 6
    {1,-m+1,-m}, // 7
    {-m,-m-1,-1}, // 8
  };
  std::vector<std::vector<int>> dir = {
    {1,2,4,8,16,32,64,128}, // 0
    {1,2,4,64,128}, // 1
    {4,8,16,32,64}, // 2
    {1,16,32,64,128}, // 3
    {1,64,128}, // 4
    {16,32,64}, // 5
    {1,2,4,8,16}, // 6
    {1,2,4}, // 7
    {4,8,16}, // 8
  };
  std::vector<int> backdir(129,0);
  backdir[1] = 16;
  backdir[2] = 32;
  backdir[4] = 64;
  backdir[8] = 128;
  backdir[16] = 1;
  backdir[32] = 2;
  backdir[64] = 4;
  backdir[128] = 8;
  std::vector<int> offset(129,0); // map from ESRI direction to index offset
  offset[1] = 1;
  offset[2] = -m+1;
  offset[4] = -m;
  offset[8] = -m-1;
  offset[16] = -1;
  offset[32] = m-1;
  offset[64] = m;
  offset[128] = m+1;

  // allocate output image initialized with 1
  RasterDataSet<int> result_image(image.originLong(),image.originLat(),image.dLong(),image.dLat(),
                                  image.sizeLong(),image.sizeLat(),1,1);
  auto& result = result_image.data();
  std::cout << "result=" << result.data() << std::endl;

  // helper data giving distance from the root
  std::vector<int> distance(N,0);
  std::cout << "distance=" << distance.data() << std::endl;

  // we need to climb up the tree to compute the distance from the outlet
  std::vector<int> front; // pixels to treat in next round
  // determine roots, i.e. pixels not draining to another pixel
  for (int index=0; index<N; index++)
    {
      if (data[index]<=0) {front.push_back(index); continue;}
      int i=index%m; // column
      int j=index/m; // row
      if (i==0 && (data[index]>=8 && data[index]<=32))
        {
          // drains to outside pixel left
          front.push_back(index);
          continue;
        }
      if (i==m-1 && (data[index]<=2 || data[index]>=128))
        {
          // drains to outside pixel right
          front.push_back(index);
          continue;
        }
      if (j==0 && (data[index]>=2 && data[index]<=8))
        {
          // drains to outside pixel below
          front.push_back(index);
          continue;
        }
      if (j==n-1 && (data[index]>=32 && data[index]<=128))
        {
          // drains to outside pixel above
          front.push_back(index);
          continue;
        }
    }

  // climb up the tree and determine distance
  int round=1;
  while (front.size()>0)
    {
      std::cout << front.size() << " pixels in front" << std::endl;
      std::vector<int> newfront; // pixels to treat in next round
      for (auto index : front)
        {
          // get the current pixel
          int i=index%m; // column
          int j=index/m; // row
          int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0)); // pixel boundary type
          
          // loop through all neighbors of current pixel
          for (int k=0; k<neighbors[c].size(); k++)
            {
              int index_neighbor=index+neighbors[c][k];
              if (data[index_neighbor]==backdir[dir[c][k]])
                {
                  // this neighbor drains to the current pixel
                  if (distance[index_neighbor]>0)
                    {
                      std::cout << "ERROR pixel was already visited distance=" << distance[index_neighbor] << std::endl;
                      exit(1);
                    }
                  distance[index_neighbor] = distance[index]+1;
		  if (distance[index_neighbor]!=round)
		  {
			  std::cout << "ERROR: distance not equal to round" << std::endl;
			  exit(1);
		  }
                  newfront.push_back(index_neighbor);
                }
            }
        }
      //std::cout << newfront.size() << " pixels in new front" << std::endl;
      front = newfront;
      round++;
    }

  // sort pixels by distance
  std::cout << "sort by distance" << std::endl;
  using X = std::pair<int,int>;
  std::vector<X> x; // pixels sorted by height 
  for (int index=0; index<N; index++)
    x.push_back(X(distance[index],index));
  std::sort(x.begin(),x.end());

  // now compute accumulation value
  std::cout << "accumulate" << std::endl;
  for (auto it=x.rbegin(); it!=x.rend(); ++it) // this is important that we proceed from high to low
    {
      int index = it->second; // retrieve index of pixel
      if (index<0 || index>=N)
      {
	      std::cout << "index out of range " << index << std::endl;
	      exit(1);
      }
      int dir = data[index]; // get drainage direction of that pixel
      int index_neighbor = index + offset[dir];
      if (distance[index]>0)
      	result[index_neighbor] += result[index];
    }
  std::cout << "finished accumulate" << std::endl;
  
  return result_image;
}


// computes image with drainage direction from an image with heigh values
template<typename T>
RasterDataSet<int> drainage_direction (const RasterDataSet<T>& image)
{
  // analyse size etc of input image
  auto& data=image.data(); // access to raw image data
  int N=data.size(); // number of pixels in image
  int m=image.sizeLong(); // pixels per row
  int n=image.sizeLat(); // pixels per column
  std::cout << "direction: N=" << N << " m=" << m << " n=" << n << std::endl;

  /* definition of drainage direction (ESRI definition)
   *
   * 32 64 128
   * 16 0 1
   * 8 4 2
   * (aold definition:)
   * 8 1 2
   * 7 0 3
   * 6 5 4
   */
  
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

  // allocate output image
  RasterDataSet<int> result_image(image.originLong(),image.originLat(),image.dLong(),image.dLat(),
                                  image.sizeLong(),image.sizeLat(),0.0,1);
  auto& result = result_image.data();
  
  // compute accumulation and drainto fields
  std::vector<int> drainto(N,-1); // a pixel drains to the lowest neighbor that is strictly lower than itself;
                                 // 0 means there is none

  for (int index=0; index<N; ++index) // this is important that we proceed from high to low
    {
      int i=index%m; // column
      int j=index/m; // row
      int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0)); // pixel boundary type
      int minheight=data[index];
      int mink=-1;
      // determine lowest lying neighbor
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
          if (neighbors[c][mink]==1) drainto[index] = 1;
          if (neighbors[c][mink]==-m+1) drainto[index] = 2;
          if (neighbors[c][mink]==-m) drainto[index] = 4;
          if (neighbors[c][mink]==-m-1) drainto[index] = 8;
          if (neighbors[c][mink]==-1) drainto[index] = 16;
          if (neighbors[c][mink]==m-1) drainto[index] = 32;
          if (neighbors[c][mink]==m) drainto[index] = 64;
          if (neighbors[c][mink]==m+1) drainto[index] = 128;
        }
    }

  result = drainto; // copy
  return result_image;
}

// catchment that drains through a given pixel
template<typename T>
RasterDataSet<int> catchment (const RasterDataSet<T>& image, int ox, int oy)
{
  // analyse size etc of input image
  auto& data=image.data(); // access to raw image data
  int N=data.size(); // number of pixels in image
  int m=image.sizeLong(); // pixels per row
  int n=image.sizeLat(); // pixels per column
  std::cout << "catchment: N=" << N << " m=" << m << " n=" << n
            << " ox=" << ox << " oy=" << oy << std::endl;

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

  // allocate output image
  // and initialize
  RasterDataSet<int> result_image(image.originLong(),image.originLat(),image.dLong(),image.dLat(),
                                    image.sizeLong(),image.sizeLat(),0,1);
  result_image(ox,oy) = 1;
  std::size_t catchment_size = 1;
  // std::cout << "added pixel " << ox << "," << oy << std::endl; 
  auto outletheight = image(ox,oy);
  auto& result = result_image.data();
  
  // sort pixels by height
  std::cout << "sort by height" << std::endl;
  using X = std::pair<int,int>;
  std::vector<X> x; // pixels sorted by height 
  for (int j=0; j<n; j++) // loop over rows
    for (int i=0; i<m; i++) // loop over columns
      {
        int index = j*m+i;
        T value = data[index];
        x.push_back(X(value,index));
      }
  std::sort(x.begin(),x.end());

  // now build catchment by proceeding from low to high
  //std::cout << "catch" << std::endl;
  for (auto it=x.begin(); it!=x.end(); ++it) // this is important that we proceed from high to low
    {
      int height = it->first;
      if (height <= outletheight) continue;
      
      // determine lowest lying neighbor
      int index = it->second;
      int i=index%m; // column
      int j=index/m; // row
      int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0)); // pixel boundary type
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

      // if we have one proliferate catchment
      if (mink>=0)
        {
          // we have a neighbor with smaller height; accumulate this pixel to the neighbor
          int index_neighbor=index+neighbors[c][mink];
          result[index] = result[index_neighbor];
          if (result[index_neighbor]==1)
            {
              catchment_size += 1;
              // std::cout << "added pixel " << i << "," << j << std::endl;
            }
        }
    }
  std::cout << "catchment size is " << catchment_size << std::endl;
  
  return result_image;
}

// number all catchments
template<typename T>
RasterDataSet<int> catchment (const RasterDataSet<T>& image)
{
  // analyse size etc of input image
  auto& data=image.data(); // access to raw image data
  int N=data.size(); // number of pixels in image
  int m=image.sizeLong(); // pixels per row
  int n=image.sizeLat(); // pixels per column
  std::cout << "catchment: N=" << N << " m=" << m << " n=" << n
            << std::endl;

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

  // allocate output image
  // and initialize
  RasterDataSet<int> result_image(image.originLong(),image.originLat(),image.dLong(),image.dLat(),
                                    image.sizeLong(),image.sizeLat(),-1,1);
  std::size_t catchment_count = 0;
  auto& result = result_image.data();
  
  // sort pixels by height
  std::cout << "sort by height" << std::endl;
  using X = std::pair<int,int>;
  std::vector<X> x; // pixels sorted by height 
  for (int j=0; j<n; j++) // loop over rows
    for (int i=0; i<m; i++) // loop over columns
      {
        int index = j*m+i;
        T value = data[index];
        x.push_back(X(value,index));
      }
  std::sort(x.begin(),x.end());

  // now build catchment by proceeding from low to high
  std::cout << "catch" << std::endl;
  for (auto it=x.begin(); it!=x.end(); ++it) // this is important that we proceed from high to low
    {
      int height = it->first;
      
      // determine lowest lying neighbor
      int index = it->second;
      int i=index%m; // column
      int j=index/m; // row
      int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0)); // pixel boundary type
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

      if (mink<0)
        {
          // we do not have a lowest neighbor; start new catchment
          if (result[index]!=-1)
            std::cout << "this pixel " << index << " has value " << result[index] << " but expected -1" << std::endl;
          result[index] = catchment_count++;
        }
      else
        {
          // we have a neighbor with smaller height; add pixel to that catchment
          int index_neighbor=index+neighbors[c][mink];
          if (result[index_neighbor]==-1)
            std::cout << "neighbor pixel " << index_neighbor << " has unexpected value -1 " << std::endl;
          result[index] = result[index_neighbor];
        }
    }

  return result_image;
}

// number all catchments
template<typename T>
void join_catchments (RasterDataSet<T>& image, const std::set<T>& source, T target)
{
  // analyse size etc of input image
  auto& data=image.data(); // access to raw image data
  auto N=data.size(); // number of pixels in image

  for (std::size_t i=0; i<N; i++)
    if (source.count(data[i])>0)
      data[i] = target;
}


// image with drainage direction
template<typename T>
RasterDataSet<float> curvature (const RasterDataSet<T>& image)
{
  // analyse size etc of input image
  auto& data=image.data(); // access to raw image data
  int N=data.size(); // number of pixels in image
  int m=image.sizeLong(); // pixels per row
  int n=image.sizeLat(); // pixels per column
  std::cout << "curvature: N=" << N << " m=" << m << " n=" << n
            << " dLong=" << image.dLong() << " dLat=" << image.dLat()
            << std::endl;

  // allocate output image
  RasterDataSet<float> result_image(image.originLong(),image.originLat(),image.dLong(),image.dLat(),
                                    image.sizeLong(),image.sizeLat(),0.0,1);
  auto& result = result_image.data();

  // compute curvature
  for (int j=1; j<n-1; j++)
    for (int i=1; i<m-1; i++)
      {
        int index = j*m+i;
        result[index] = data[index-m]+data[index+m]+data[index-1]+data[index+1]-4.0*data[index];
      }

  // return image
  return result_image;
}

// image with drainage direction
template<typename T>
RasterDataSet<float> derivative_x (const RasterDataSet<T>& image)
{
  // analyse size etc of input image
  auto& data=image.data(); // access to raw image data
  int N=data.size(); // number of pixels in image
  int m=image.sizeLong(); // pixels per row
  int n=image.sizeLat(); // pixels per column
  std::cout << "derivative_x: N=" << N << " m=" << m << " n=" << n
            << " dLong=" << image.dLong() << " dLat=" << image.dLat()
            << std::endl;

  // allocate output image
  RasterDataSet<float> result_image(image.originLong(),image.originLat(),image.dLong(),image.dLat(),
                                    image.sizeLong(),image.sizeLat(),0.0,1);
  auto& result = result_image.data();

  // compute curvature
  for (int j=1; j<n-1; j++)
    for (int i=1; i<m-1; i++)
      {
        int index = j*m+i;
        result[index] = 0.5*(data[index+1]-data[index-1]);
      }

  // return image
  return result_image;
}

// image with drainage direction
template<typename T>
RasterDataSet<float> derivative_y (const RasterDataSet<T>& image)
{
  // analyse size etc of input image
  auto& data=image.data(); // access to raw image data
  int N=data.size(); // number of pixels in image
  int m=image.sizeLong(); // pixels per row
  int n=image.sizeLat(); // pixels per column
  std::cout << "derivative_x: N=" << N << " m=" << m << " n=" << n
            << " dLong=" << image.dLong() << " dLat=" << image.dLat()
            << std::endl;

  // allocate output image
  RasterDataSet<float> result_image(image.originLong(),image.originLat(),image.dLong(),image.dLat(),
                                    image.sizeLong(),image.sizeLat(),0.0,1);
  auto& result = result_image.data();

  // compute curvature
  for (int j=1; j<n-1; j++)
    for (int i=1; i<m-1; i++)
      {
        int index = j*m+i;
        result[index] = 0.5*(data[index+m]-data[index-m]);
      }

  // return image
  return result_image;
}

// image with drainage direction
template<typename T>
RasterDataSet<float> flatness (const RasterDataSet<T>& image)
{
  std::cout << "flatness" << std::endl;
  
  // analyse size etc of input image
  auto& data=image.data(); // access to raw image data
  int N=data.size(); // number of pixels in image
  int m=image.sizeLong(); // pixels per row
  int n=image.sizeLat(); // pixels per column
  std::cout << "derivative_x: N=" << N << " m=" << m << " n=" << n
            << " dLong=" << image.dLong() << " dLat=" << image.dLat()
            << std::endl;

  // allocate output image
  RasterDataSet<float> result_image(image.originLong(),image.originLat(),image.dLong(),image.dLat(),
                                    image.sizeLong(),image.sizeLat(),0.0,1);
  auto& result = result_image.data();

  // compute curvature
  for (int j=1; j<n-1; j++)
    for (int i=1; i<m-1; i++)
      {
        int index = j*m+i;
        result[index] = 0.5*std::max(abs(data[index+m]-data[index-m]),abs(data[index+1]-data[index-1]));
      }

  // return image
  return result_image;
}

// determine pixels with groundwater reservoir
// 0) initialize pixels with -1
// 1) the pixels with accumulation >= minaccumulation get value 0
// 2) if a pixel neighbors a pixel with value>=0 and |gradx|+|grady|<=maxslope it gets value of neighbor+1
template<typename T1, typename T2>
RasterDataSet<int> groundwaterbody (const RasterDataSet<T1>& _accumulation,
                                const RasterDataSet<T2>& _flatness,
                                float minaccumulation,
                                float maxslope
                                )
{
  // analyse size etc of input image
  auto& accumulation=_accumulation.data(); // access to raw image data
  auto& flatness=_flatness.data(); // access to raw image data
  int N=accumulation.size(); // number of pixels in image
  int m=_accumulation.sizeLong(); // pixels per row
  int n=_accumulation.sizeLat(); // pixels per column
  std::cout << "groundwater: N=" << N << " m=" << m << " n=" << n << std::endl;

  // parameters
  
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

  // allocate output image
  RasterDataSet<int> result_image(_accumulation.originLong(),_accumulation.originLat(),_accumulation.dLong(),_accumulation.dLat(),
                                  _accumulation.sizeLong(),_accumulation.sizeLat(),-1,1); // fill with -1
  auto& result = result_image.data();
  
  // initialize
  std::vector<int> x; // pixels to be looked at in next round
  int d=0;     // distance of pixels in this round
  for (int j=0; j<n; j++) // loop over rows
    for (int i=0; i<m; i++) // loop over columns
      {
        int index = j*m+i;
        if (accumulation[index]>=minaccumulation)
          {
            result[index] = 0; // a valley center pixel
            int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0)); // pixel boundary type
            for (int k=0; k<neighbors[c].size(); k++) 
              {
                // look at all neighbors
                int index_neighbor=index+neighbors[c][k];
                if (result[index_neighbor]<0 && accumulation[index_neighbor]<minaccumulation)
                  x.push_back(index_neighbor);
              }
          }
      }
  // now pixels are either 0 when in the valley center or -1
  std::cout << "distance=0 initialized" << std::endl;

  // extend to flat regions connected to valley center
  while (x.size()>0)
    {
      d += 1; // distance in this round
      std::cout << "distance=" << d << " looking at " << x.size() << " pixels" << std::endl; 
      std::vector<int> y; // determine pixels to consider in next round
      for (auto index : x) // look at all candidates
        if (result[index]<0 && flatness[index]<=maxslope)
          {
            result[index]=d; // do not consider that pixel anymore
            int i=index%m; // column
            int j=index/m; // row
            int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0)); // pixel boundary type
            for (int k=0; k<neighbors[c].size(); k++)
              {
                int index_neighbor=index+neighbors[c][k];
                if (result[index_neighbor]<0)
                  y.push_back(index_neighbor);
              }
          }
      x = y; // go to next iteration
    }
  // now pixels are
  // -1 : no groundwater
  // 0 : valley center
  // d> distance from valley center


  // produce another field with
  // distance from boundary of groundwater
  int dmax = 1000000000;
  std::vector<int> distance(N,dmax);

  // find pixels which are in groundwater body but have a non-groundwater neighbor
  x.clear();
  for (int index=0; index<N; index++)
    {
      if (result[index]<0) continue; // then this is not a groundwater pixel
      
      // get the current pixel
      int i=index%m; // column
      int j=index/m; // row
      int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0)); // pixel boundary type

      // inspect neighbors
      for (int k=0; k<neighbors[c].size(); k++)
        {
          int index_neighbor=index+neighbors[c][k];
          if (result[index_neighbor]<0)
            {
              distance[index] = 0; // this pixel is at the boundary of the valley
              x.push_back(index); // we should consider its neighbors in the next round
              break; // we have treated that pixel
            }
        }
    }
  // distance has either dmax or 0
  // x contains a list of pixels which are at the valley boundary (have value 0)

  d = 0;
  while (x.size()>0)
    {
      d += 1; // distance in this round
      std::cout << "distance=" << d << " looking at " << x.size() << " pixels" << std::endl; 
      std::vector<int> y; // determine pixels to consider in next round
      for (auto index : x) // look at all candidates
        {
          // this index has been treated already; we need to look at it neighbors
          int i=index%m; // column
          int j=index/m; // row
          int c = 3*(2*(j+1==n)+(j==0))+(2*(i+1==m)+(i==0)); // pixel boundary type
          for (int k=0; k<neighbors[c].size(); k++)
            {
              int index_neighbor=index+neighbors[c][k];
              if (result[index_neighbor]<0) continue; // not groundwater, skip it
              if (distance[index_neighbor]>distance[index]+1)
                {
                  distance[index_neighbor] = distance[index]+1;
                  y.push_back(index_neighbor);
                }
            }
        }
      x = y; // go to next iteration
    }
  // distance of a pixel is either dmax or it is 0<=d<dmax and dmax is the distance from the boundary

    // fill in pixels at the boundary that connect gw pixels over edges
  for (int j=1; j<n-1; j++) // loop over inner rows
    for (int i=1; i<m-1; i++) // loop over inner columns
      {
        int index = j*m+i;
        if (result[index]<0) 
          continue; // then this is not a groundwater pixel
        // count connections
        int connections=0;
        if (result[index-1]>=0) connections++; // connects to left
        if (result[index+1]>=0) connections++; // connects to right
        if (result[index-m]>=0) connections++; // connects to down
        if (result[index+m]>=0) connections++; // connects to up
        if (connections>=2) continue; // is already connected
        // pixel is groundwater but is connected to less than two 
        // other groundwater pixels via edges, so connect it
        // find candidate neighbors
        std::set<int> s;
        if (result[index+m]<0 && (result[index+m-1]>=0 || result[index+m+1]>=0))
          s.insert(index+m);
        if (result[index-m]<0 && (result[index-m-1]>=0 || result[index-m+1]>=0))
          s.insert(index-m);
        if (result[index-1]<0 && (result[index+m-1]>=0 || result[index-m-1]>=0))
          s.insert(index-1);
        if (result[index+1]<0 && (result[index+m+1]>=0 || result[index-m+1]>=0))
          s.insert(index+1);
        // now select the flatest ones
        int maxflatness = 10000;
        int bestindex = -1;
        for (auto ix : s) // look at all candidates
          if (flatness[ix]<maxflatness)
            {
              maxflatness = flatness[ix];
              bestindex = ix;
            }
        if (bestindex<0) continue; // ups, nothing to connect
        result[bestindex] = result[index]+1;
        distance[bestindex] = -1;
        s.erase(bestindex);
        connections++;
        if (connections>=2) continue;
        // one more time
        maxflatness = 10000;
        bestindex = -1;
        for (auto ix : s) // look at all candidates
          if (flatness[ix]<maxflatness)
            {
              maxflatness = flatness[ix];
              bestindex = ix;
            }
        if (bestindex<0) continue; // ups, nothing to connect
        result[bestindex] = result[index]+1;
        distance[bestindex] = -1;
      }

  // produce final result
  for (int index=0; index<N; index++)
    {
      auto value=distance[index];
      if (value==dmax)
        distance[index] = -1; // no groundwater
      else if (value==-1)
        distance[index] = 0; // additional fill pixels
      else distance[index] = value+1; // real groundwater pixels
    }

  result = distance;
  return result_image;
}

#endif
