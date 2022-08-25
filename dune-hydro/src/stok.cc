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
#include <fstream> 
#include <sstream>
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
// dune-istl includes
#include<dune/istl/matrixmarket.hh>
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
#include<dune/pdelab/gridfunctionspace/subspace.hh>
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
#include <dune/vtk/pvdwriter.hh>

// dune-nonlinopt
#include<dune/nonlinopt/nonlinopt.hh>

// include stuff from dune-hydro
#include<dune/hydro/nonlineardiffusionfv.hh>
#include<dune/hydro/nonlineardiffusionfv_velocity.hh>
#include<dune/hydro/geotiffreader.hh>
#include<dune/hydro/netcdfreader.hh>
#include<dune/hydro/rasterdataset.hh>
#include<dune/hydro/rasterdatasetalgorithms.hh>
#include<dune/hydro/nonlineardiffusioncoupledfv.hh>

// check if a file exists
bool exists_file (const std::string& name)
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

/** \brief get number of day that time value t in seconds lies in
 *
 * It is defined as follows t is in day n if
 *
 * (n==0 && t <= 86400+epsilon) || (n*86400+epsilon < t and t <= (n+1)*86400+epsilon) 
 *
 * use epsilon = 1e-4
 */
template<typename T>
int day_from_time (T t)
{
  T epsilon = 1e-4;
  if (t-epsilon<=86400.0) return 0;
  return std::floor((t-epsilon)/86400.0);
}


//********************************************************************************
/** \brief A class holding and comparing measurement data

 */
//********************************************************************************

class DataAnalysis
{
public:

  // load data file
  std::map<int,double> load_measurement_file (std::string filename)
  {
    std::map<int,double> measurement;
    std::ifstream myfile(filename);
    if (myfile.is_open())
      {
        while ( !myfile.eof() )
          {
            std::string date;
            myfile >> date;
            if (myfile.eof()) break;
            double day;
            myfile >> day;
            if (myfile.eof()) break;
            double discharge;
            myfile >> discharge;
            if (myfile.eof()) break;
            double accumulated_discharge;
            myfile >> accumulated_discharge;
            measurement[std::round(day-1.0)] = discharge;
            //std::cout << "DATA " << date << " " << day*86400.0 << " " << discharge << std::endl;
          }
        myfile.close();
        std::cout << "load_measurement_file: " << measurement.size() << " items read" << std::endl;
      }
    return measurement;
  }

  std::vector<std::pair<double,double>> load_data_file (std::string filename)
  {
    std::vector<std::pair<double,double>> data;
    std::ifstream myfile(filename);
    if (myfile.is_open())
      {
        while ( !myfile.eof() )
          {
            std::string s;
            myfile >> s;
            if (myfile.eof()) break;
            double time;
            myfile >> time;
            time *= 86400.0;
            if (myfile.eof()) break;
            double discharge_s;
            myfile >> discharge_s;
            if (myfile.eof()) break;
            double discharge_g;
            myfile >> discharge_g;
            if (myfile.eof()) break;
            double accumulated_discharge_s;
            myfile >> accumulated_discharge_s;
            if (myfile.eof()) break;
            double accumulated_discharge_g;
            myfile >> accumulated_discharge_g;
            if (myfile.eof()) break;
            double catchment_s;
            myfile >> catchment_s;
            if (myfile.eof()) break;
            double catchment_g;
            myfile >> catchment_g;
            if (myfile.eof()) break;

            data.push_back(std::pair<double,double>(time,discharge_s));
            //std::cout << "DATA " << date << " " << day*86400.0 << " " << discharge << std::endl;
          }
        myfile.close();
        std::cout << "load_data_file: " << data.size() << " items read" << std::endl;
      }
    return data;
  }

  // compute distance of daily averages
  double l2_distance (const std::map<int,double>& measurement, const std::map<int,double>& simulation, int from=0)
  {
    // simulation and measurements are supposed to be one value per day!
    double error = 0.0;
    
    // the simulation results <t_i,d_i> means constant discharge d_i on the interval (t_{i-1},t_i)
    for (auto xs : simulation)
      {
        auto iterm = measurement.find(xs.first);
        if (iterm!=measurement.end())
          error += (iterm->second-xs.second)*(iterm->second-xs.second);
      }
    std::cout << "l2_distance: " << simulation.size() << " items processed" << std::endl;
    return std::sqrt(error);
  }

  // compute dayly averages and return them in a map
  template<typename C>
  std::map<int,double> daily_averages (const C& simulation)
  {
    std::map<int,double> averages;
    std::map<int,double> times;
    auto iter1 = simulation.begin();
    if (iter1==simulation.end()) return averages;
    auto iter2 = iter1;
    iter2++;
    int n=0;
    while (iter2!=simulation.end())
      {
        auto time = iter2->first;
        auto dt = time - iter1->first;
        auto discharge = iter2->second;
        int day = 1+day_from_time(time);
        auto iter = averages.find(day);
        if (iter==averages.end())
          {
            // new entry
            averages[day] = discharge*dt;
            times[day] = dt;
          }
        else
          {
            // new entry
            averages[day] += discharge*dt;
            times[day] += dt;
          }
        iter1 = iter2;
        iter2++;
        n++;
      }
    for (auto iter=averages.begin(); iter!=averages.end(); ++iter)
      iter->second /= times[iter->first];
    std::cout << "daily_averages: " << n << " items processed and " << averages.size() << " days generated" << std::endl;
    return averages;
  }
};

//********************************************************************************
/** \brief Model class provides parameters to a shallow scalar and coupled flow model 

   \tparam Number : a number type
 */
//********************************************************************************

template<typename Number>
class Model
{
  enum {dim=2};
  Dune::ParameterTree ptree;
  double dx,dy,ox,oy;
  Dune::FieldVector<Number,dim> L;
  Dune::FieldVector<Number,dim> H;
  std::array<int,dim> N;
  Number precipitationrate;
  Number initialheight;
  RasterDataSet<short> elevation_raster;
  RasterDataSet<float> gwdepth_raster;
  RasterDataSet<short> flow_direction;
  RasterDataSet<float> flatnessmap;
  RasterDataSet<int> accumulationmap;
  RasterDataSet<int> groundwaterdepthmap;
  Number eps,A,B;
  double time;
  RasterDataSet<short> ga;
  std::vector<RasterDataSet<short>> sc;
  std::vector<RasterDataSet<float>> temp;
  std::vector<RasterDataSet<float>> precip;
  std::vector<RasterDataSet<float>> lst;
  Number k_surface_value, porosity, conductivity, L_i, L_e, C, ds, d0,dmax;
  std::vector<std::size_t> days_per_month;
  std::string path_to_data;

  void add_data (int year, int days)
  {
    // temperature
    std::string combinedname = path_to_data + "DISTRIBUTED_T_STATION/";
    std::cout << "=====> READING TEMPERATURE DATA for " << year << std::endl;
    int n=0;
    for (int month=1; month<=12; month++)
      for (int day=1; day<=days_per_month[month-1]; day++)
        if (n<days)
          {
            char filename[512];
            sprintf(filename,"%04d-%02d-%02d.tif",year,month,day);
            GeoTIFFImage<double> image(filename,combinedname,1,2);
            temp.push_back(RasterDataSet<float>(ox+0.5*dx,oy+0.5*dy,dx,dy,N[0],N[1],-9.9999,1));
            temp.back().paste(image);
            n += 1;
          }
    std::cout << temp.size() << " temperature files read successfully" << std::endl;

    // precipitation
    n = 0;
    std::cout << "=====> READING PRECIPITATION DATA for " << year << std::endl;
    combinedname = path_to_data + "PRECP_DISTRIBUTED_2018-19/";
    for (int month=1; month<=12; month++)
      for (int day=1; day<=days_per_month[month-1]; day++)
        if (n<days)
          {
            char filename[512];
            sprintf(filename,"%04d-%02d-%02d.tif",year,month,day);
            GeoTIFFImage<double> image(filename,combinedname,1,2);
            precip.push_back(RasterDataSet<float>(ox+0.5*dx,oy+0.5*dy,dx,dy,N[0],N[1],0.0,1));
            precip.back().paste(image);
            n += 1;
          }
    std::cout << temp.size() << " precipitation files read successfully" << std::endl;

    // land surface temperature
    std::cout << "=====> READING LAND SURFACE TEMPERATURE DATA for " << year << std::endl;
    combinedname = path_to_data + "LST_ALL_STOK_2018_2019/";
    n = 0;
    for (int day=1; day<=365; day++)
      if (n<days)
        {
          char filename[512];
          sprintf(filename,"FILLED_MODIS_Regression_LST_USING_AIRT_STOK_%04d_%d.tif",year,day);
          bool success=true;
          try {
            GeoTIFFImage<double> image(filename,combinedname,1,2);
            lst.push_back(RasterDataSet<float>(ox+0.5*dx,oy+0.5*dy,dx,dy,N[0],N[1],-9.9999,1));
            lst.back().paste(image);
          }
          catch (GeoTIFFReaderException &e){
            success = false; // could not read file, try other data type
          }
          if (!success)
            {
              GeoTIFFImage<float> image(filename,combinedname,1,2);
              lst.push_back(RasterDataSet<float>(ox+0.5*dx,oy+0.5*dy,dx,dy,N[0],N[1],-9.9999,1));
              lst.back().paste(image);
            }
          n += 1;
        }
    std::cout << temp.size() << " land surface temperature files read successfully" << std::endl;
  }
  
public:
  //! Constructor gets parameter tree to read ini file parameters
  Model (Dune::ParameterTree ptree_)
    : ptree(ptree_), eps(1e-2)
  {
    std::vector<std::size_t> dpm = {31,28,31,30,31,30,31,31,30,31,30,31};
    days_per_month = dpm;

    // path to all the data
    path_to_data = ptree.get<std::string>("problem.path_to_data");

    // read hydrosheds elevation data
    GeoTIFFImage<GInt16> n30e075("n30e075_con.tif",path_to_data,1,2);
    dx=std::abs(n30e075.dLong());
    dy=std::abs(n30e075.dLat());
    ox = ptree.get("problem.ox",(double)77.3987); 
    oy = ptree.get("problem.oy",(double)33.89998); 
    N[0] = ptree.get("problem.NX",(int)320);
    N[1] = ptree.get("problem.NY",(int)320);
    H[0] = 90.0;
    H[1] = 90.0;
    L[0] = N[0]*H[0];
    L[1] = N[1]*H[1];
    elevation_raster = RasterDataSet<short>(ox+0.5*dx,oy+0.5*dy,dx,dy,N[0],N[1],0,1);
    elevation_raster.paste(n30e075);
    elevation_raster = outlier(elevation_raster,2,2000);

    // modify pixel values
    elevation_raster(130,151) = 4029;
    std::cout << "MODIFY PIXEL " << 130 << "," << 151 << " from " << 4023 <<  " to " << 4029 << std::endl;
    elevation_raster.data()[33056] = 4850; 
    std::cout << "MODIFY PIXEL " << 33056%N[0] << "," << 33056/N[0] << " from " << 4857 <<  " to " << 4850 << std::endl;
    elevation_raster.data()[44280] = 4206;
    std::cout << "MODIFY PIXEL " << 44280%N[0] << "," << 44280/N[0] << " from " << 4209 <<  " to " << 4206 << std::endl;
    elevation_raster.data()[44921] = 4189;
    std::cout << "MODIFY PIXEL " << 44921%N[0] << "," << 44921/N[0] << " from " << 4191 <<  " to " << 4189 << std::endl;
    elevation_raster.data()[45884] = 4150;
    std::cout << "MODIFY PIXEL " << 45884%N[0] << "," << 45884/N[0] << " from " << 4153 <<  " to " << 4150 << std::endl;
    elevation_raster.data()[48769] = 4022;
    std::cout << "MODIFY PIXEL " << 48769%N[0] << "," << 48769/N[0] << " from " << 4028 <<  " to " << 4022 << std::endl;
    elevation_raster.data()[49409] = 4018;
    std::cout << "MODIFY PIXEL " << 49409%N[0] << "," << 49409/N[0] << " from " << 4026 <<  " to " << 4018 << std::endl;
    elevation_raster.data()[50048] = 4000;
    std::cout << "MODIFY PIXEL " << 50048%N[0] << "," << 50048/N[0] << " from " << 4002 <<  " to " << 4000 << std::endl;
    elevation_raster.data()[51012] = 3945;
    std::cout << "MODIFY PIXEL " << 51012%N[0] << "," << 51012/N[0] << " from " << 3947 <<  " to " << 3945 << std::endl;
    elevation_raster.data()[50693] = 3927;
    std::cout << "MODIFY PIXEL " << 50693%N[0] << "," << 50693/N[0] << " from " << 3930 <<  " to " << 3927 << std::endl;
    elevation_raster.data()[52309] = 3778;
    std::cout << "MODIFY PIXEL " << 52309%N[0] << "," << 52309/N[0] << " from " << 3781 <<  " to " << 3778 << std::endl;
    elevation_raster.data()[53914] = 3729;
    std::cout << "MODIFY PIXEL " << 53914%N[0] << "," << 53914/N[0] << " from " << 3731 <<  " to " << 3729 << std::endl;
    elevation_raster.data()[55517] = 3693;
    std::cout << "MODIFY PIXEL " << 55517%N[0] << "," << 55517/N[0] << " from " << 3698 <<  " to " << 3693 << std::endl;
    elevation_raster.data()[35968] = 4481;
    std::cout << "MODIFY PIXEL " << 35968%N[0] << "," << 35968/N[0] << " from " << 4490 <<  " to " << 4481 << std::endl;

    // other parameters
    k_surface_value = ptree.get<double>("problem.k_surface");
    porosity = ptree.get<double>("problem.porosity");
    conductivity = ptree.get<double>("problem.conductivity");
    L_i = ptree.get<double>("problem.L_i");
    L_e = ptree.get<double>("problem.L_e");
    C = ptree.get<double>("problem.C");
    ds = ptree.get<double>("problem.ds");
    d0 = ptree.get<double>("problem.d0");
    dmax = ptree.get<double>("problem.dmax");
    std::cout << "MODEL PARAMS k_surface=" << k_surface_value << " d0=" << d0 << " L_i=" << L_i << std::endl;

    // flow direction
    GeoTIFFImage<GInt16> n30e075dir("n30e075_dir.tif",path_to_data,1,2);
    flow_direction = RasterDataSet<short>(ox+0.5*dx,oy+0.5*dy,dx,dy,N[0],N[1],0,1);
    flow_direction.paste(n30e075dir);

    // glacierized area
    std::string combinedname = path_to_data + "GLACERISED_AREA_RASTER/";
    GeoTIFFImage<GByte> gaimage("STOK_GLAC_2019_FINAL90m_A.tif",combinedname,1,2);
    ga = RasterDataSet<short>(ox+0.5*dx,oy+0.5*dy,dx,dy,N[0],N[1],0,1);
    ga.paste(gaimage);
    
    // snow cover data
    combinedname = path_to_data + "MODIS_SCA_2018_2019_DAILY_90m/";
    for (int year=2018; year<=2019; year++)
      for (int i=1; i<=365; i++)
        {
          char filename[512];
          sprintf(filename,"MOYD10A1GL06_Maximum_Snow_Extent_%04d%03d.tif",year,i);
          GeoTIFFImage<GByte> image(filename,combinedname,1,2);
          sc.push_back(RasterDataSet<short>(ox+0.5*dx,oy+0.5*dy,dx,dy,N[0],N[1],0,1));
          sc.back().paste(image);
        }

    // groundwater depth
    accumulationmap = accumulation_from_direction(flow_direction);
    flatnessmap = flatness(elevation_raster);
    groundwaterdepthmap = groundwaterbody(accumulationmap,flatnessmap,100,22);
    gwdepth_raster = RasterDataSet<float>(ox+0.5*dx,oy+0.5*dy,dx,dy,N[0],N[1],0.0,1);
    for (std::size_t i=0; i<gwdepth_raster.data().size(); i++)
      {
        if (groundwaterdepthmap.data()[i]==-1)
          gwdepth_raster.data()[i] = 0.25; // bare rock pixel
        else if (groundwaterdepthmap.data()[i]==0)
          gwdepth_raster.data()[i] = 1.0; // the fill pixels
        else
          gwdepth_raster.data()[i] = std::min(ds*groundwaterdepthmap.data()[i]+d0,dmax);
        if (ga.data()[i]==200)
          gwdepth_raster.data()[i] = 0.5; // under the glacier pixel
      }

    add_data(2018,365);
    add_data(2019,273);

    // compute total melt
    double total_melt=0.0; // in cubicmeter
    for (std::size_t d=0; d<temp.size(); d++) // loop over all days
      for (std::size_t i=0; i<temp[0].data().size(); i++)
        {
          auto snowcover_value = sc[d].data()[i];
          auto temperature_value = temp[d].data()[i];
          if (snowcover_value==100 && temperature_value>0.0)
            total_melt += temperature_value*3.9*1e-3*90*90;
        }
    std::cout << "total melt in " << temp.size() << " days is " << total_melt << " cubicmeter" << std::endl;

    // read precipitation rate
    precipitationrate = ptree.get("problem.precipitationrate",(Number)0.0);
    initialheight = ptree.get("problem.initialheight",(Number)0.0);

    // regularization parameter
    A = 2.5/(2.0*std::sqrt(eps));
    B = 0.5/(2.0*pow(eps,2.5));
    std::cout << "REGULARIZATION gradient eps=" << eps << " A=" << A << " B=" << B << std::endl;
  }

  //*****************************************************************************
  // export data
  //*****************************************************************************

  const Dune::FieldVector<Number,dim>& length () const
  {
    return L;
  }

  const std::array<int,dim>& cells () const
  {
    return N;
  }

  const RasterDataSet<short>& elevation ()
  {
    return elevation_raster;
  }

  const RasterDataSet<float>& depth ()
  {
    return gwdepth_raster;
  }

  const RasterDataSet<int>& gwdepth ()
  {
    return groundwaterdepthmap;
  }

  const RasterDataSet<short>& direction ()
  {
    return flow_direction;
  }

  const RasterDataSet<short>& glacier ()
  {
    return ga;
  }

  const RasterDataSet<short>& snowcover (std::size_t i)
  {
    i = std::min(i,sc.size()-1);
    return sc[i];
  }

  const RasterDataSet<float>& temperature (std::size_t i)
  {
    i = std::min(i,temp.size()-1);
    return temp[i];
  }

  const RasterDataSet<float>& precipitation (std::size_t i)
  {
    i = std::min(i,precip.size()-1);
    return precip[i];
  }

  const RasterDataSet<float>& land_surface_temperature (std::size_t i)
  {
    i = std::min(i,precip.size()-1);
    return lst[i];
  }
  
  //*****************************************************************************
  // scalar equation
  //*****************************************************************************

  //! storage term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number storage (const E& e, const X& x) const
  {
    return 1.0;
  }
  
  //! bottom position in global coordinates (function b above)
  template<typename E, typename X>
  /**
     \param e      : codim 0 entity
     \param xlocal : local coordinate
   */
  Number bathymmetry (const E& e, const X& xlocal) const
  {
    auto x = e.geometry().center();
    int cellx = std::floor(x[0]/H[0]);
    int celly = std::floor(x[1]/H[1]);
    return elevation_raster(cellx,celly);
  }

  //! nonlinearity to be evaluated at a point on a face seperating two elements
  /**
     \param b_inside     : bottom position on inside element
     \param b_outside    : bottom position on outside element
     \param u_inside     : surface position on inside element
     \param u_outside    : surface position on outside element
     \param vn           : velocity in normal direction on face 
   */
  template<typename U, typename VN>
  Number phi (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    // upwind evaluation of height
    Number u_upwind = (vn>=0) ? u_inside : u_outside;
    
    // evaluation of (u-b)
    Number bmax=std::max(b_inside,b_outside);
    Number h_upwind = std::max(u_upwind-bmax,0.0);
    
    // regularization of gradient and evaluation
    Number gradterm =  (std::abs(vn)>=eps) ? 1.0/std::sqrt(std::abs(vn)) : (A-B*vn*vn);

    return h_upwind*gradterm;
  }

  template<typename U, typename VN>
  Number dphi_du_inside (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    // check if we depend on u_inside
    if (vn<0) return 0.0;

    // upwind evaluation of height
    Number u_upwind = u_inside;

    // evaluation of (u-b)
    if (u_upwind-std::max(b_inside,b_outside)<0.0) return 0.0;
    
    // regularization of gradient and evaluation
    Number gradterm =  (std::abs(vn)>=eps) ? 1.0/std::sqrt(std::abs(vn)) : (A-B*vn*vn);

    return gradterm;
  }

  template<typename U, typename VN>
  Number dphi_du_outside (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    // check if we depend on u_inside
    if (vn>=0) return 0.0;

    // upwind evaluation of height
    Number u_upwind = u_outside;

    // evaluation of (u-b)
    if (u_upwind-std::max(b_inside,b_outside)<0.0) return 0.0;
    
    // regularization of gradient and evaluation
    Number gradterm =  (std::abs(vn)>=eps) ? 1.0/std::sqrt(std::abs(vn)) : (A-B*vn*vn);

    return gradterm;
  }

  template<typename U, typename VN>
  Number dphi_dvn (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    // upwind evaluation of height
    Number u_upwind = (vn>=0) ? u_inside : u_outside;

    // evaluation of (u-b)
    Number bmax=std::max(b_inside,b_outside);
    Number h_upwind = std::max(u_upwind-bmax,0.0);
    
    // regularization of gradient and evaluation
    Number sign = (vn>=0) ? 1.0 : -1.0;
    if (std::abs(vn)>=eps)
      return -h_upwind*0.5/(std::abs(vn)*std::sqrt(std::abs(vn)))*sign;
    else
      return -h_upwind*B*2.0*vn;
  }
  
  //! permeability value at position in a cell
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number k (const E& e, const X& x) const
  {
    return k_surface_value;
  }
  
  //! source/sink term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number f (const E& e, const X& xlocal) const
  {
    auto x = e.geometry().center();
    double mm_per_day_to_meter_per_second = 1e-3/86400.0;

    int cellx = std::floor(x[0]/H[0]);
    int celly = std::floor(x[1]/H[1]);
    int timestep = day_from_time(time);
    timestep = std::min(timestep,(int)sc.size()-1);
    timestep = std::min(timestep,(int)temp.size()-1);
    
    auto snowcover_value = sc[timestep](cellx,celly);
    auto temperature_value = temp[timestep](cellx,celly);
    //auto lst_value = lst[timestep](cellx,celly);
    auto precipitation_value = precip[timestep](cellx,celly) * mm_per_day_to_meter_per_second;
    auto glacier_value = ga(cellx,celly);

    Number result=0.0;

    // snow melt
    if (snowcover_value==100 && temperature_value>0.0) // if there is snow and T>0
      {
        result += temperature_value*3.1*mm_per_day_to_meter_per_second;
      }

    // ice melt
    if (glacier_value==200 && temperature_value>0.0) // if there is ice and T>0
      {
        result += temperature_value*5.9*mm_per_day_to_meter_per_second;
      }

    // additional snow melt induced by rain
    if (snowcover_value==100 && temperature_value>=1.0) // if there is snow and T>=1
      {
        result += 4.2 * precipitation_value * temperature_value / (335.0*0.96);
      }
    
    // additional ice melt induced by rain
    if (glacier_value==200 && temperature_value>=1.0) // if there is ice and T>=1
      {
        result +=  4.2 * precipitation_value * temperature_value / 335.0;
      }

    // rain that does not fall on ice or snow
    if (snowcover_value!=100 && glacier_value!=200 && temperature_value>=1.0 && precipitation_value>0.0) // if there is no ice and no snow
      {
        result +=  precipitation_value;
      }

    return result;
  }

  //! boundary condition type function (true = Dirichlet)
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  bool b (const I& i, const X& x) const
  {
    return true;
  }

  //! Dirichlet extension
  /** Provides Dirichlet boundary values as well as initial value!
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number g (const E& e, const X& xlocal) const
  {
    auto height = bathymmetry(e,xlocal);
    auto x = e.geometry().global(xlocal);
    Number eps=1e-4;
    if (x[0]>eps && x[0]<L[0]-eps && x[1]>eps && x[1]<L[1]-eps)
      return height+initialheight;
    else
      return height;
  }

  //! Neumann boundary condition
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  Number j (const I& i, const X& x) const
  {
    return 0.0;
  }

  //*****************************************************************************
  // surface equation
  // delegate to the functions above
  //*****************************************************************************

  //! storage term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number storage_surface (const E& e, const X& x) const
  {
    return storage(e,x);
  }

  //! bottom position in global coordinates (function b above)
  /**
     \param e      : codim 0 entity
     \param xlocal : local coordinate
   */
  template<typename E, typename X>
  Number bathymmetry_surface (const E& e, const X& xlocal) const
  {
    return bathymmetry(e,xlocal);
  }

  //! nonlinearity to be evaluated at a point on a face seperating two elements
  /**
     \param b_inside     : bottom position on inside element
     \param b_outside    : bottom position on outside element
     \param u_inside     : surface position on inside element
     \param u_outside    : surface position on outside element
     \param vn           : velocity in normal direction on face 
   */
  template<typename U, typename VN>
  Number phi_surface (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    return phi(b_inside,b_outside,u_inside,u_outside,vn);
  }

  //! permeability value at position in a cell
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number k_surface (const E& e, const X& x) const
  {
    return k(e,x);
  }
  
  //! source/sink term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number f_surface (const E& e, const X& x) const
  {
    return f(e,x);
  }
  
  //! boundary condition type function (true = Dirichlet)
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  bool b_surface (const I& i, const X& x) const
  {
    return b(i,x);
  }

  //! Dirichlet extension
  /** Provides Dirichlet boundary values as well as initial value!
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number g_surface (const E& e, const X& x) const
  {
    return g(e,x);
  }

  //! Neumann boundary condition
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  Number j_surface (const I& i, const X& x) const
  {
    return j(i,x);
  } 

  //*****************************************************************************
  // groundwater equation
  //*****************************************************************************

  //! storage term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number storage_groundwater (const E& e, const X& x) const
  {
    return porosity;
  }

  //! bottom position in global coordinates (function b above)
  /**
     \param e      : codim 0 entity
     \param xlocal : local coordinate
   */
  template<typename E, typename X>
  Number bathymmetry_groundwater (const E& e, const X& xlocal) const
  {
    auto x = e.geometry().center();
    int cellx = std::floor(x[0]/H[0]);
    int celly = std::floor(x[1]/H[1]);
    return elevation_raster(cellx,celly)-gwdepth_raster(cellx,celly);
  }

  //! nonlinearity to be evaluated at a point on a face seperating two elements
  /**
     \param b_inside     : bottom position on inside element
     \param b_outside    : bottom position on outside element
     \param u_inside     : surface position on inside element
     \param u_outside    : surface position on outside element
     \param vn           : velocity in normal direction on face 
   */
  template<typename U, typename VN>
  Number phi_groundwater (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    // upwind evaluation of height
    Number u_upwind = (vn>=0) ? u_inside : u_outside;
    
    // evaluation of (u-b)
    Number bmax=std::max(b_inside,b_outside);
    Number h_upwind = std::max(u_upwind-bmax,0.0);
    
    return h_upwind;
  }

  //! permeability value at position in a cell
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number k_groundwater (const E& e, const X& xlocal) const
  {
    auto x = e.geometry().center();
    int cellx = std::floor(x[0]/H[0]);
    int celly = std::floor(x[1]/H[1]);
    if (gwdepth_raster(cellx,celly)<1.0)
      return 0.0;
    else
      return conductivity;
  }
  
  //! source/sink term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number f_groundwater (const E& e, const X& x) const
  {
    return 0.0;
  }
  
  //! boundary condition type function (true = Dirichlet)
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  bool b_groundwater (const I& i, const X& x) const
  {
    return true;
  }

  //! Dirichlet extension
  /** Provides Dirichlet boundary values as well as initial value!
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number g_groundwater (const E& e, const X& x) const
  {
    return bathymmetry_groundwater(e,x);
  }

  //! Neumann boundary condition
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  Number j_groundwater (const I& i, const X& x) const
  {
    return 0.0;
  } 

  //*****************************************************************************
  // coupling
  //*****************************************************************************

  // Exchange term
  template<typename E, typename X>
  Number q (const E& e, const X& xlocal, Number u_surface, Number u_groundwater, Number bath_surf) const
  {
    auto x = e.geometry().center();
    int cellx = std::floor(x[0]/H[0]);
    int celly = std::floor(x[1]/H[1]);

    if (gwdepth_raster(cellx,celly)<1.0)
      return 0.0; // no exchange with rock or under the glacier

    Number value1 = smoothmax0(u_surface-bath_surf,0.01); // std::max(u_surface-bath_surf,0.0);
    Number value2 = smoothmax0(u_surface-u_groundwater,0.01); // std::max(u_surface-u_groundwater,0.0);
    Number value3 = smoothmax0(u_groundwater-u_surface,0.01); // std::max(u_groundwater - u_surface,0.0);
    Number firstterm  = L_i*(value1/(C + value1))*value2;
    Number secondterm = L_e*value3;

    return firstterm-secondterm;
  }

  //*****************************************************************************
  // time setter
  //*****************************************************************************
  
  //! set time for subsequent evaluation
  /**
     \param t value of time for subsequent evaluations
   */
  void setTime (double t)
  {
    time = t;
  }
};


//********************************************************************************
// scalar driver function solving the problem
//********************************************************************************

void driverFVNewNewton (Dune::ParameterTree& ptree)
{
  // dimension and important types
  const int dim = 2;
  using RF = double;               // type for computations

  using Grid = Dune::YaspGrid<dim>;
  using DF = Grid::ctype;
  using GV = Grid::LeafGridView;

  using FEM = Dune::PDELab::P0LocalFiniteElementMap<RF,RF,dim>;
  using CON = Dune::PDELab::P0ParallelConstraints;
  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;

  using Z = Dune::PDELab::Backend::Vector<GFS,RF>;

  using ZDGF = Dune::PDELab::DiscreteGridFunction<GFS,Z>;
  using VTKF = Dune::PDELab::VTKGridFunctionAdapter<ZDGF>;

  using MODEL = Model<RF>;

  // make time
  RF time = ptree.get("problem.tstart",(RF)0.0);

  // make model
  MODEL model(ptree);
  model.setTime(time);

  // make YaspGrid
  int overlap = ptree.get("problem.overlap",(int)1);
  std::shared_ptr<Grid> gridp = std::make_shared<Grid>(model.length(),model.cells(),std::bitset<dim>(0ULL),overlap);
  GV gv=gridp->leafGridView();

  // Make grid function space
  FEM fem(Dune::GeometryTypes::cube(dim));
  CON con;
  GFS gfs(gv,fem,con);
  gfs.name("Vh");

  // A coefficient vector
  Z z(gfs); // initial value

  // Make a grid function out of it
  ZDGF zdgf(gfs,z);

  // make user functions and set initial time
  auto glambda = [&](const auto& e, const auto& x)
    {return model.g(e,x);};
  auto g = Dune::PDELab::
    makeInstationaryGridFunctionFromCallable(gv,glambda,model);
  auto blambda = [&](const auto& i, const auto& x)
    {return model.b(i,x);};
  auto b = Dune::PDELab::
    makeBoundaryConditionFromCallable(gv,blambda);
  auto bathymmetrylambda = [&](const auto& e, const auto& x)
    {return model.bathymmetry(e,x);};
  auto bathymmetrygf = Dune::PDELab::
    makeInstationaryGridFunctionFromCallable(gv,bathymmetrylambda,model);
  Z bathymmetry(gfs);
  Dune::PDELab::interpolate(bathymmetrygf,gfs,bathymmetry);
  ZDGF bdgf(gfs,bathymmetry);

  // Assemble constraints
  typedef typename GFS::template
    ConstraintsContainer<RF>::Type CC;
  CC cc; cc.clear();
  Dune::PDELab::constraints(b,gfs,cc); // assemble constraints
  std::cout << "constrained dofs=" << cc.size() << " of "
            << gfs.globalSize() << std::endl;

  // initialize simulation time,  the coefficient vector
  Dune::PDELab::interpolate(g,gfs,z);

  // Make instationary grid operator
  typedef NonlinearDiffusionFV<MODEL,GFS> LOP;
  LOP lop(model,gfs);
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(2*dim+1); // guess nonzeros per row
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,CC,CC> GO0;
  GO0 go0(gfs,cc,gfs,cc,lop,mbe);

  typedef FVL2<MODEL> TLOP;
  TLOP tlop(model);
  typedef Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,RF,RF,RF,CC,CC> GO1;
  GO1 go1(gfs,cc,gfs,cc,tlop,mbe);

  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
  IGO igo(go0,go1);
  igo.divideMassTermByDeltaT();
  //igo.multiplySpatialTermByDeltaT();

  // Select a linear solver backend
  //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<IGO> LS;
  //typedef Dune::PDELab::ISTLBackend_SEQ_UMFPack LS;
  //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILU0 LS;
  //LS ls(1000,0);
  typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS;
  LS ls (gfs,100,1);
  auto params = ls.parameters();
  params.setCoarsenTarget(60000);
  ls.setParameters(params);
  
  // using LS = Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC>;
  // LS ls(gfs,cc,500,5,2);

  // solve nonlinear problem
  typedef Dune::PDELab::NewtonMethod<IGO,LS> PDESOLVER;
  PDESOLVER pdesolver(igo,ls);
  pdesolver.setParameters(ptree.sub("newton"));
  if (gv.comm().rank()!=0)
    pdesolver.setVerbosityLevel(0);

  pdesolver.setReassembleThreshold(0.0);
  pdesolver.setVerbosityLevel(2);
  pdesolver.setAbsoluteLimit(1e-7);
  pdesolver.setReduction(1e-6);
  pdesolver.setMinLinearReduction(1e-4);
  pdesolver.setUseMaxNorm(true);
  //pdesolver.setMaxIterations(10);
  //pdesolver.setLineSearchMaxIterations(8);
  //pdesolver.setLineSearchStrategy("hackbuschReuskenAcceptBest");
  //pdesolver.setLineSearchStrategy("hackbuschReusken");

  // select and prepare time-stepping scheme
  Dune::PDELab::OneStepThetaParameter<RF> method1(1.0);
  Dune::PDELab::Alexander2Parameter<RF> method2;
  Dune::PDELab::Alexander3Parameter<RF> method3;
  int torder = ptree.get("method.torder",(int)1);
  Dune::PDELab::TimeSteppingParameterInterface<RF>*
    pmethod=&method1;
  if (torder==1) pmethod = &method1;
  if (torder==2) pmethod = &method2;
  if (torder==3) pmethod = &method3;
  if (torder<1||torder>3)
    std::cout<<"torder not in [1,3]"<<std::endl;
  typedef Dune::PDELab::OneStepMethod<RF,IGO,PDESOLVER,Z,Z> OSM;
  OSM  osm(*pmethod,igo,pdesolver);
  osm.setVerbosityLevel(2);

  // velocity postprocessing
  Z znew(z);
  ZDGF znewdgf(gfs,znew);
  using VeloDGF = NonlinearDiffusionFVVelocity<MODEL,GFS,Z>;
  using VeloVTKF = Dune::PDELab::VTKGridFunctionAdapter<VeloDGF>;
  VeloDGF velodgf(model,gfs,znew,bathymmetry);

  // residuals
  Z rspace(gfs); rspace = 0.0;  ZDGF rspacedgf(gfs,rspace);
  Z rtime(gfs); rtime = 0.0;    ZDGF rtimedgf(gfs,rtime);
  Z rtotal(gfs); rtotal = 0.0;  ZDGF rtotaldgf(gfs,rtotal);
  Z zdiff(gfs); zdiff = 0.0;    ZDGF zdiffdgf(gfs,zdiff);

  // prepare VTK writer and write first file
  std::string filename=ptree.get("output.filename","output");
  // struct stat st;
  // if( stat( filename.c_str(), &st ) != 0 )
  //   {
  //     int stat = 0;
  //     stat = mkdir(filename.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
  //     if( stat != 0 && stat != -1)
  //       std::cout << "Error: Cannot create directory "
  //                 << filename << std::endl;
  //   }

  // data analysis and visualization of bathymmetry
  // auto catchment_image = catchment(model.elevation(),161,179);
  auto catchment_image = catchment(model.elevation());
  std::set<int> joiners;
  joiners.insert(1343);
  joiners.insert(1325);
  joiners.insert(1341);
  joiners.insert(1340);
  joiners.insert(1370);
  joiners.insert(1377);
  joiners.insert(1368);
  //joiners.insert(1362);
  join_catchments(catchment_image,joiners,513);
  Z catchment_vec(gfs);
  ZDGF catchment_dgf(gfs,catchment_vec);
  for (size_t i=0; i<catchment_image.data().size(); i++)
    Dune::PDELab::Backend::native(catchment_vec)[i][0] = catchment_image.data()[i];

  // accumulation computed from elevation
  auto accumulation_image1 = drainage_accumulation(model.elevation());
  Z accumulation_vec1(gfs);
  ZDGF accumulation_dgf1(gfs,accumulation_vec1);
  for (size_t i=0; i<accumulation_image1.data().size(); i++)
    Dune::PDELab::Backend::native(accumulation_vec1)[i][0] = accumulation_image1.data()[i];

  // accumulation computed from hydroSHEDS directions
  auto accumulation_image2 = accumulation_from_direction(model.direction());
  Z accumulation_vec2(gfs);
  ZDGF accumulation_dgf2(gfs,accumulation_vec2);
  for (size_t i=0; i<accumulation_image2.data().size(); i++)
    Dune::PDELab::Backend::native(accumulation_vec2)[i][0] = accumulation_image2.data()[i];
  
  // curvature
  auto curvature_image = curvature(model.elevation());
  Z curvature_vec(gfs);
  ZDGF curvature_dgf(gfs,curvature_vec);
  for (size_t i=0; i<curvature_image.data().size(); i++)
    Dune::PDELab::Backend::native(curvature_vec)[i][0] = curvature_image.data()[i];

  // gradient
  auto flatness_image = flatness(model.elevation());
  Z flatness_vec(gfs);
  ZDGF flatness_dgf(gfs,flatness_vec);
  for (size_t i=0; i<flatness_image.data().size(); i++)
    Dune::PDELab::Backend::native(flatness_vec)[i][0] = flatness_image.data()[i];

  // groundwaterbody
  auto gwb_image = model.gwdepth();
  Z gwb_vec(gfs);
  ZDGF gwb_dgf(gfs,gwb_vec);
  for (size_t i=0; i<gwb_image.data().size(); i++)
    Dune::PDELab::Backend::native(gwb_vec)[i][0] = gwb_image.data()[i];

  // depth
  auto depth_image = model.depth();
  Z depth_vec(gfs);
  ZDGF depth_dgf(gfs,depth_vec);
  for (size_t i=0; i<depth_image.data().size(); i++)
    Dune::PDELab::Backend::native(depth_vec)[i][0] = depth_image.data()[i];

  // snow cover
  Z snowcover_vec(gfs);
  ZDGF snowcover_dgf(gfs,snowcover_vec);
  for (size_t i=0; i<model.snowcover(0).data().size(); i++)
    Dune::PDELab::Backend::native(snowcover_vec)[i][0] = model.snowcover(0).data()[i];

  // glacier cover
  Z glacier_vec(gfs);
  ZDGF glacier_dgf(gfs,glacier_vec);
  for (size_t i=0; i<model.glacier().data().size(); i++)
    Dune::PDELab::Backend::native(glacier_vec)[i][0] = model.glacier().data()[i];

  // temperature
  Z temperature_vec(gfs);
  ZDGF temperature_dgf(gfs,temperature_vec);
  for (size_t i=0; i<model.temperature(0).data().size(); i++)
    Dune::PDELab::Backend::native(temperature_vec)[i][0] = model.temperature(0).data()[i];

  // precipitation
  Z precipitation_vec(gfs);
  ZDGF precipitation_dgf(gfs,precipitation_vec);
  for (size_t i=0; i<model.precipitation(0).data().size(); i++)
    Dune::PDELab::Backend::native(precipitation_vec)[i][0] = model.precipitation(0).data()[i];

  // lst
  Z lst_vec(gfs);
  ZDGF lst_dgf(gfs,lst_vec);
  for (size_t i=0; i<model.land_surface_temperature(0).data().size(); i++)
    Dune::PDELab::Backend::native(lst_vec)[i][0] = model.land_surface_temperature(0).data()[i];

  // source term
  Z source_vec(gfs);
  ZDGF source_dgf(gfs,source_vec);
  auto& indexset = gv.indexSet();
  Dune::FieldVector<double,2> x0(0.5);
  for (const auto& e : elements(gv))
    Dune::PDELab::Backend::native(source_vec)[indexset.index(e)][0] = model.f(e,x0);

  // cell index
  Z index_vec(gfs);
  ZDGF index_dgf(gfs,index_vec);
  for (size_t i=0; i< Dune::PDELab::Backend::native(index_vec).size(); i++)
    Dune::PDELab::Backend::native(index_vec)[i][0] = i;

  // compute amount of water in catchment
  double water_in_catchment = 0.0;
  for (size_t i=0; i<catchment_image.data().size(); i++)
    if (catchment_image.data()[i]==513)
      {
        double height = Dune::PDELab::Backend::native(znew)[i][0] - Dune::PDELab::Backend::native(bathymmetry)[i][0];
        water_in_catchment += height*90.0*90.0;
      }
  std::cout << "BALANCE " << time << " " << water_in_catchment << std::endl;
  
  // new writer
  using Writer = Dune::VtkImageDataWriter<GV>;
  Dune::PvdWriter<Writer> pvdWriter(gv, Dune::Vtk::FormatTypes::COMPRESSED, Dune::Vtk::DataTypes::FLOAT32);
  pvdWriter.addCellData(std::make_shared<VTKF>(bdgf,"bathymmetry"));
  pvdWriter.addCellData(std::make_shared<VTKF>(index_dgf,"dune index"));
  pvdWriter.addCellData(std::make_shared<VTKF>(catchment_dgf,"catchment"));
  pvdWriter.addCellData(std::make_shared<VTKF>(snowcover_dgf,"snow cover"));
  pvdWriter.addCellData(std::make_shared<VTKF>(glacier_dgf,"glacier"));
  pvdWriter.addCellData(std::make_shared<VTKF>(temperature_dgf,"temperature"));
  pvdWriter.addCellData(std::make_shared<VTKF>(lst_dgf,"land surface temperature"));
  pvdWriter.addCellData(std::make_shared<VTKF>(precipitation_dgf,"precipitation"));
  pvdWriter.addCellData(std::make_shared<VTKF>(source_dgf,"source"));
  pvdWriter.addCellData(std::make_shared<VTKF>(accumulation_dgf1,"accumulation1"));
  pvdWriter.addCellData(std::make_shared<VTKF>(accumulation_dgf2,"accumulation2"));
  pvdWriter.addCellData(std::make_shared<VTKF>(gwb_dgf,"groundwater"));
  pvdWriter.addCellData(std::make_shared<VTKF>(depth_dgf,"depth groundwater"));
  pvdWriter.addCellData(std::make_shared<VTKF>(curvature_dgf,"curvature"));
  pvdWriter.addCellData(std::make_shared<VTKF>(flatness_dgf,"flatness"));
  pvdWriter.addCellData(std::make_shared<VTKF>(znewdgf,"solution"));
  pvdWriter.addCellData(std::make_shared<VeloVTKF>(velodgf,"velocity"));
  std::string fullfilename = filename + ".vti";
  pvdWriter.writeTimestep(time,fullfilename);
  
  // time loop
  RF T = ptree.get("problem.T",(RF)1.0);
  RF dt = ptree.get("method.dt",(RF)0.1);
  RF timestepmax =  ptree.get("method.dtmax",dt);
  int every = ptree.get("output.every",(int)1);
  int step=0;
  int increased_step = 0;
  int decreased_step = 0;
  double eps=1e-6;
  auto factor = 1.0/std::sqrt(2.0);
  double accumulated_discharge = 0.0;
  double accumulated_source_term = 0.0;
  int next_day_boundary = 1;
  
  while (time<T-1e-8)
    {
      // adjust time step to hit day boundary
      // what is the next date boundary to hit?
      double time_to_hit = next_day_boundary*86400.0;

      if (time+dt>time_to_hit-1e-5)
        {
          dt = time_to_hit-time;
        }
      else if (time+2*dt>time_to_hit-1e-5)
        {
          dt = 0.5*(time_to_hit-time);
        }        
      // assemble constraints for new time step
      model.setTime(time+dt);

      // do time step
      try {
        znew = z;
        osm.apply(time,dt,z,znew);
      }
      catch (Dune::PDELab::LineSearchError) {
        dt *= factor;
        if (dt<eps) throw;
        decreased_step = step;
        if (gv.comm().rank()==0)
          std::cout << "STEP CONTROL: " << " decreased in step " << step+1 << " newdt: " << dt << std::endl;
        continue;
      }
      catch (Dune::PDELab::TerminateError) {
        dt *= factor;
        if (dt<eps) throw;
        decreased_step = step;
        if (gv.comm().rank()==0)
          std::cout << "STEP CONTROL: " << " decreased in step " << step+1 << " newdt: " << dt << std::endl;
        continue;
      }
      catch (Dune::PDELab::NewtonLinearSolverError) {
        dt *= factor;
        if (dt<eps) throw;
        decreased_step = step;
        if (gv.comm().rank()==0)
          std::cout << "STEP CONTROL: " << " decreased in step " << step+1 << " newdt: " << dt << std::endl;
        continue;
      }
      catch (Dune::ISTLError) {
        dt *= factor;
        if (dt<eps) throw;
        decreased_step = step;
        if (gv.comm().rank()==0)
          std::cout << "STEP CONTROL: " << " decreased in step " << step+1 << " newdt: " << dt << std::endl;
        continue;
      }

      // analyze residual
      rspace = 0.0;
      go0.residual(znew,rspace);
      Dune::PDELab::set_constrained_dofs(cc,0.0,rspace);
      rtime=0.0;
      zdiff=znew;
      zdiff -= z;
      go1.residual(zdiff,rtime);
      Dune::PDELab::set_constrained_dofs(cc,0.0,rtime);
      //rtime *= 1.0/dt; // this seems to be there already!
      rtotal = rtime;
      rtotal += rspace;

      auto norm2_time = std::sqrt(gv.comm().sum(rtime.two_norm2()));
      auto norm2_space = std::sqrt(gv.comm().sum(rspace.two_norm2()));
      auto norm2_total = std::sqrt(gv.comm().sum(rtotal.two_norm2()));
      if (gv.comm().rank()==0)
        {
          std::cout << "||rtime ||_2=" << norm2_time  << std::endl;
          std::cout << "||rspace||_2=" << norm2_space  << std::endl;
          std::cout << "||rtotal||_2=" << norm2_total  << std::endl;
        }
      auto norm8_time = gv.comm().max(rtime.infinity_norm());
      auto norm8_space = gv.comm().max(rspace.infinity_norm());
      auto norm8_total = gv.comm().max(rtotal.infinity_norm());
      if (gv.comm().rank()==0)
        {
          std::cout << "||rtime ||_8=" << norm8_time  << std::endl;
          std::cout << "||rspace||_8=" << norm8_space  << std::endl;
          std::cout << "||rtotal||_8=" << norm8_total  << std::endl;
        }

      // analyze change in time step
      double max_change = 0.0;
      auto& indexset = gv.indexSet();
      for (const auto& e : elements(gv,Dune::Partitions::interior)) {
        auto i = indexset.index(e);
        auto bi = Dune::PDELab::Backend::native(bathymmetry)[i][0];
        auto ziold = Dune::PDELab::Backend::native(z)[i][0];
        auto zinew = Dune::PDELab::Backend::native(znew)[i][0];
        max_change = std::max(max_change,std::abs(ziold-zinew));
      }
      gv.comm().max(max_change);
      if (gv.comm().rank()==0)
        std::cout << "MAXCHANGE " << max_change << " " << time+dt << " " << dt << std::endl;
        
      // accept time step
      velodgf.update();
      time+=dt;
      step++;
      z = znew;

      // update temperature and snow cover
      int timestep = day_from_time(time);
      timestep = std::min(timestep,637);
      for (size_t i=0; i<model.snowcover(0).data().size(); i++)
        Dune::PDELab::Backend::native(snowcover_vec)[i][0] = model.snowcover(timestep).data()[i];
      for (size_t i=0; i<model.temperature(0).data().size(); i++)
        Dune::PDELab::Backend::native(temperature_vec)[i][0] = model.temperature(timestep).data()[i];
      for (size_t i=0; i<model.precipitation(0).data().size(); i++)
        Dune::PDELab::Backend::native(precipitation_vec)[i][0] = model.precipitation(timestep).data()[i];
      for (size_t i=0; i<model.land_surface_temperature(0).data().size(); i++)
        Dune::PDELab::Backend::native(lst_vec)[i][0] = model.land_surface_temperature(timestep).data()[i];

      // compute amount of water in catchment
      double water_in_catchment = 0.0;
      for (size_t i=0; i<catchment_image.data().size(); i++)
        if (catchment_image.data()[i]==513)
          {
            double height = Dune::PDELab::Backend::native(znew)[i][0] - Dune::PDELab::Backend::native(bathymmetry)[i][0];
            water_in_catchment += height*90.0*90.0;
          }

      // update source term
      model.setTime(time);
      if (time/86400.0>365.0 && time/86400.0<375.0)
        accumulated_source_term = water_in_catchment;
      for (const auto& e : elements(gv))
        {
          Dune::PDELab::Backend::native(source_vec)[indexset.index(e)][0] = model.f(e,x0);
          accumulated_source_term += dt*90*90*Dune::PDELab::Backend::native(source_vec)[indexset.index(e)][0];
        }

      std::cout << "BALANCE " << time
                << " " << water_in_catchment
                << " " << accumulated_source_term
                << std::endl;

      // probe for discharge output
      using VeloProbe = Dune::PDELab::GridFunctionProbe<VeloDGF>;
      using Probe = Dune::PDELab::GridFunctionProbe<ZDGF>;
      using Domain = typename VeloDGF::Traits::DomainType;
      using VeloRange = typename VeloDGF::Traits::RangeType;
      using Range = typename ZDGF::Traits::RangeType;
      Domain point;
      point[0] = 14550.0; point[1] = 16190.0; // (161,179)
      VeloProbe veloprobe(velodgf,point);
      Probe solprobe(znewdgf,point);
      Probe bathyprobe(bdgf,point);
      VeloRange veloresult; veloprobe.eval_all(veloresult);
      Range solresult; solprobe.eval_all(solresult);
      Range bathyresult; bathyprobe.eval_all(bathyresult);
      // accumulated_discharge += dt*(solresult[0] - bathyresult[0])*90.0*veloresult[1];
      accumulated_discharge += dt*90.0*veloresult[1];
      if (time/86400.0>365.0 && time/86400.0<375.0) accumulated_discharge = 0.0;
      if (gv.comm().rank()==0)
        std::cout << "DISCHARGE "
                  << time/86400.0
                  << " " << veloresult[0]
                  << " " << veloresult[1]
                  << " " << solresult[0] - bathyresult[0]
                  << " " << 90.0*veloresult[1]
                  << " " << accumulated_discharge
                  << " " << water_in_catchment 
                  << " " << accumulated_source_term
                  << " " << (accumulated_source_term - (water_in_catchment + accumulated_discharge))/(accumulated_source_term+1e-6)*100.0
                  << std::endl;

      // determine maximum velicity value in any edge
      auto maxv = maxvelocity(velodgf);
      if (gv.comm().rank()==0)
        std::cout << "MAXVELOCITY " << maxv << std::endl;
      
      // output to VTK file
      if (step%every==0)
        {
          pvdWriter.writeTimestep(time,fullfilename);
          // vtkSequenceWriter.write(time,Dune::VTK::appendedraw);
          if (gv.comm().rank()==0)
            std::cout << "WRITING VTK OUTPUT " << step << " " << time << " " << dt << std::endl;
        }

      // increase timestep
      if (step-decreased_step>5 && step-increased_step>5)
        {
          double newdt = std::min(timestepmax,std::sqrt(2.0)*dt);
          if (newdt>dt)
            {
              increased_step = step;
              if (gv.comm().rank()==0)
                std::cout << "STEP CONTROL: " << " increased in step " << step << " newdt: " << newdt << std::endl;
              dt = newdt;
            }
        }
    }
}

namespace Dune {
  namespace PDELab {

   /**
     * @brief Overlapping parallel CG solver with UMFPack preconditioner
     * @tparam GFS The Type of the GridFunctionSpace.
     * @tparam CC The Type of the Constraints Container.
     */
    template<class GFS, class CC>
    class ISTLBackend_OVLP_BCGS_UMFPack
      : public ISTLBackend_OVLP_UMFPack_Base<GFS,CC,Dune::BiCGSTABSolver>
    {
    public:

      /*! \brief make a linear solver object

        \param[in] gfs_ a grid function space
        \param[in] cc_ a constraints object
        \param[in] maxiter_ maximum number of iterations to do
        \param[in] verbose_ print messages if true
      */
      ISTLBackend_OVLP_BCGS_UMFPack (const GFS& gfs_, const CC& cc_,
                                              unsigned maxiter_=5000,
                                              int verbose_=1)
        : ISTLBackend_OVLP_UMFPack_Base<GFS,CC,Dune::BiCGSTABSolver>(gfs_,cc_,maxiter_,verbose_)
      {}
    };

  }
}

//********************************************************************************
// coupled driver function solving the problem
//********************************************************************************

std::map<int,double> driverCoupledFVNewNewton (Dune::ParameterTree& ptree)
{
  // dimension and important types
  const int dim = 2;
  using RF = double;               // type for computations

  using Grid = Dune::YaspGrid<dim>;
  using DF = Grid::ctype;
  using GV = Grid::LeafGridView;

  using FEM = Dune::PDELab::P0LocalFiniteElementMap<RF,RF,dim>;
  using CON = Dune::PDELab::P0ParallelConstraints;
  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
  using OrderingTag = Dune::PDELab::EntityBlockedOrderingTag;
  using GFS2 = Dune::PDELab::PowerGridFunctionSpace<GFS,2,VBE,OrderingTag>;
  using Path0 = Dune::TypeTree::HybridTreePath<Dune::index_constant<0>>;
  using Path1 = Dune::TypeTree::HybridTreePath<Dune::index_constant<1>>;
  using GFS20 = Dune::PDELab::GridFunctionSubSpace<GFS2,Path0>;
  using GFS21 = Dune::PDELab::GridFunctionSubSpace<GFS2,Path1>;
  using CC2 = GFS2::template ConstraintsContainer<RF>::Type;

  using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  using Z2 = Dune::PDELab::Backend::Vector<GFS2,RF>;

  using ZDGF = Dune::PDELab::DiscreteGridFunction<GFS,Z>;
  using VTKF = Dune::PDELab::VTKGridFunctionAdapter<ZDGF>;

  using Z2DGF = Dune::PDELab::DiscreteGridFunction<GFS2,Z2>;
  using Z2DGF0 = Dune::PDELab::DiscreteGridFunction<GFS20,Z2>;
  using Z2DGF1 = Dune::PDELab::DiscreteGridFunction<GFS21,Z2>;
  using VTKF20 = Dune::PDELab::VTKGridFunctionAdapter<Z2DGF0>;
  using VTKF21 = Dune::PDELab::VTKGridFunctionAdapter<Z2DGF1>;

  using MODEL = Model<RF>;

  using VeloDGFSurface = NonlinearDiffusionFVVelocitySurface<MODEL,GFS,Z>;
  using VeloVTKFSurface = Dune::PDELab::VTKGridFunctionAdapter<VeloDGFSurface>;
  using VeloDGFGroundwater = NonlinearDiffusionFVVelocityGroundwater<MODEL,GFS,Z>;
  using VeloVTKFGroundwater = Dune::PDELab::VTKGridFunctionAdapter<VeloDGFGroundwater>;

  using namespace Dune::Indices;

  // read time to start with; we could start with 2019 immediately
  RF time = ptree.get("problem.tstart",(RF)0.0);

  // make model
  MODEL model(ptree);
  model.setTime(time);

  // make YaspGrid
  std::cout << "make grid" << std::endl;
  int overlap = ptree.get("problem.overlap",(int)1);
  std::shared_ptr<Grid> gridp = std::make_shared<Grid>(model.length(),model.cells(),std::bitset<dim>(0ULL),overlap);
  GV gv=gridp->leafGridView();

  // Make grid function space
  std::cout << "make gfs" << std::endl;
  FEM fem(Dune::GeometryTypes::cube(dim));
  CON con;
  GFS gfs(gv,fem,con);
  gfs.name("Vh");

  // make combine GFS with two components
  GFS2 gfs2(gfs);
  gfs2.name("VCoupled");
  gfs2.child(_0).name("surface");
  gfs2.child(_1).name("groundwater");
  GFS20 gfs20(gfs2); // subspace for surface component
  GFS21 gfs21(gfs2); // subspace for subsurface component

  // A coefficient vector
  Z2 z2(gfs2);   // value at old time step
  Z2 z2new(z2); // value at new time step

  // Make discrete grid functions for components
  Z2DGF0 z2newdgf0(gfs20,z2new);
  Z2DGF1 z2newdgf1(gfs21,z2new);

  // bathymmetries
  std::cout << "make bathymmetries" << std::endl;
  auto bathymmetrylambda_surface = [&](const auto& e, const auto& x)
    {return model.bathymmetry_surface(e,x);};
  auto bathymmetrysurfacegf = Dune::PDELab::
    makeInstationaryGridFunctionFromCallable(gv,bathymmetrylambda_surface,model);
  Z bathymmetry_surface(gfs);
  Dune::PDELab::interpolate(bathymmetrysurfacegf,gfs,bathymmetry_surface);
  ZDGF bsurfacedgf(gfs,bathymmetry_surface);
  auto bathymmetrylambda_groundwater = [&](const auto& e, const auto& x)
    {return model.bathymmetry_groundwater(e,x);};
  auto bathymmetrygroundwatergf = Dune::PDELab::
    makeInstationaryGridFunctionFromCallable(gv,bathymmetrylambda_groundwater,model);
  Z bathymmetry_groundwater(gfs);
  Dune::PDELab::interpolate(bathymmetrygroundwatergf,gfs,bathymmetry_groundwater);
  ZDGF bgroundwaterdgf(gfs,bathymmetry_groundwater);

  // initial condition
  std::cout << "make initial conditions" << std::endl;
  auto glambda_surface = [&](const auto& e, const auto& x)
    {return model.g_surface(e,x);};
  auto g_surface = Dune::PDELab::
    makeInstationaryGridFunctionFromCallable(gv,glambda_surface,model);
  auto glambda_groundwater = [&](const auto& e, const auto& x)
    {return model.g_groundwater(e,x);};
  auto g_groundwater = Dune::PDELab::
    makeInstationaryGridFunctionFromCallable(gv,glambda_groundwater,model);

  // read initial value from file if true else start with initial value in class
  if (ptree.get<bool>("output.loadinitial"))
    {
      std::string loadfilename=ptree.get<std::string>("output.loadfilename");
      auto s = loadfilename + "_" + std::to_string(gv.comm().rank()) + ".mm";
      loadMatrixMarket(Dune::PDELab::Backend::native(z2),s);
    }
  else
    {
      Dune::PDELab::interpolate(g_surface,gfs20,z2);
      Dune::PDELab::interpolate(g_groundwater,gfs21,z2);
    }
  z2new = z2;
  
  // boundary condition type
  std::cout << "make boundary conditions" << std::endl;
  auto blambda_surface = [&](const auto& i, const auto& x)
    {return model.b_surface(i,x);};
  auto b_surface = Dune::PDELab::
    makeBoundaryConditionFromCallable(gv,blambda_surface);
  auto blambda_groundwater = [&](const auto& i, const auto& x)
    {return model.b_groundwater(i,x);};
  auto b_groundwater = Dune::PDELab::
    makeBoundaryConditionFromCallable(gv,blambda_groundwater);
  using B = Dune::PDELab::CompositeConstraintsParameters<decltype(b_surface),decltype(b_groundwater)>;
  B b(b_surface, b_groundwater);
  
  // Assemble constraints
  std::cout << "make constraints" << std::endl;
  CC2 cc2; cc2.clear();
  Dune::PDELab::constraints(b,gfs2,cc2); // assemble constraints
  std::cout << "constrained dofs=" << cc2.size() << " of " << gfs2.globalSize() << std::endl;

  // Make instationary grid operator
  std::cout << "make operators" << std::endl;
  typedef NonlinearDiffusionCoupledFV<MODEL,GFS2> LOP;
  LOP lop(model,gfs2);
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(2*dim+1); // guess nonzeros per row
  typedef Dune::PDELab::GridOperator<GFS2,GFS2,LOP,MBE,RF,RF,RF,CC2,CC2> GO0;
  GO0 go0(gfs2,cc2,gfs2,cc2,lop,mbe);

  typedef FVTemporalCoupled<MODEL,GFS2> TLOP;
  TLOP tlop(model,gfs2);
  typedef Dune::PDELab::GridOperator<GFS2,GFS2,TLOP,MBE,RF,RF,RF,CC2,CC2> GO1;
  GO1 go1(gfs2,cc2,gfs2,cc2,tlop,mbe);

  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
  IGO igo(go0,go1);
  igo.divideMassTermByDeltaT();
  //igo.multiplySpatialTermByDeltaT();

  // Select a linear solver backend
  std::cout << "make linear solver" << std::endl;
  //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<IGO> LS;
  //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILU0 LS;

  // typedef Dune::PDELab::ISTLBackend_SEQ_UMFPack LS;
  // LS ls(1000,0);

  //typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS;
  //LS ls (gfs,100,1);
  //auto params = ls.parameters();
  //params.setCoarsenTarget(60000);
  //ls.setParameters(params);
  // using LS = Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC>;
  // LS ls(gfs,cc,500,5,2);
  using LS = Dune::PDELab::ISTLBackend_OVLP_BCGS_UMFPack<GFS2,CC2>;
  LS ls(gfs2,cc2,500,2);

  // solve nonlinear problem
  std::cout << "make newton" << std::endl;
  typedef Dune::PDELab::NewtonMethod<IGO,LS> PDESOLVER;
  PDESOLVER pdesolver(igo,ls);
  pdesolver.setParameters(ptree.sub("newton"));
  if (gv.comm().rank()!=0) pdesolver.setVerbosityLevel(0);
  pdesolver.setReassembleThreshold(0.0);
  pdesolver.setVerbosityLevel(2);
  pdesolver.setAbsoluteLimit(1e-7);
  pdesolver.setReduction(1e-6);
  pdesolver.setMinLinearReduction(1e-4);
  pdesolver.setUseMaxNorm(true); 
  //pdesolver.setMaxIterations(10);
  //pdesolver.setLineSearchMaxIterations(8);
  //pdesolver.setLineSearchStrategy("hackbuschReuskenAcceptBest");
  //pdesolver.setLineSearchStrategy("hackbuschReusken");

  // select and prepare time-stepping scheme
  std::cout << "make timestepper" << std::endl;
  Dune::PDELab::OneStepThetaParameter<RF> method1(1.0);
  Dune::PDELab::Alexander2Parameter<RF> method2;
  Dune::PDELab::Alexander3Parameter<RF> method3;
  int torder = ptree.get("method.torder",(int)1);
  Dune::PDELab::TimeSteppingParameterInterface<RF>*
    pmethod=&method1;
  if (torder==1) pmethod = &method1;
  if (torder==2) pmethod = &method2;
  if (torder==3) pmethod = &method3;
  if (torder<1||torder>3)
    std::cout<<"torder not in [1,3]"<<std::endl;
  using OSM = Dune::PDELab::OneStepMethod<RF,IGO,PDESOLVER,Z2,Z2>;
  OSM  osm(*pmethod,igo,pdesolver);
  osm.setVerbosityLevel(2);

  // velocity postprocessing
  Z z_surface(gfs);
  ZDGF zdgf_surface(gfs,z_surface);
  Dune::PDELab::interpolate(z2newdgf0,gfs,z_surface);
  VeloDGFSurface velodgfsurface(model,gfs,z_surface,bathymmetry_surface);
  Z z_groundwater(gfs);
  ZDGF zdgf_groundwater(gfs,z_groundwater);
  Dune::PDELab::interpolate(z2newdgf1,gfs,z_groundwater);
  VeloDGFGroundwater velodgfgroundwater(model,gfs,z_groundwater,bathymmetry_groundwater);

  // height postprocessing
  Z h_surface(z_surface);
  ZDGF hdgf_surface(gfs,h_surface);
  h_surface -= bathymmetry_surface;
  Z h_groundwater(z_groundwater);
  ZDGF hdgf_groundwater(gfs,h_groundwater);
  h_groundwater -= bathymmetry_groundwater;

  // prepare VTK writer and write first file
  std::cout << "make data vectors" << std::endl;
  std::string filename=ptree.get("output.filename","output");
  auto& indexset = gv.indexSet();

  // data analysis and visualization of bathymmetry
  // auto catchment_image = catchment(model.elevation(),161,179);
  std::cout << "call catchment" << std::endl;
  auto catchment_image = catchment(model.elevation());
  std::set<int> joiners;
  joiners.insert(1343);
  joiners.insert(1325);
  joiners.insert(1341);
  joiners.insert(1340);
  joiners.insert(1370);
  joiners.insert(1377);
  joiners.insert(1368);
  //joiners.insert(1362);
  std::cout << "call join_catchments" << std::endl;
  join_catchments(catchment_image,joiners,513);
  Z catchment_vec(gfs);
  ZDGF catchment_dgf(gfs,catchment_vec);
  for (const auto& e : elements(gv))
    {
      auto x = e.geometry().center();
      int cellx = std::floor(x[0]/90.0);
      int celly = std::floor(x[1]/90.0);
      Dune::PDELab::Backend::native(catchment_vec)[indexset.index(e)][0] = catchment_image(cellx,celly);
    }

  // accumulation computed from hydroSHEDS directions
  std::cout << "call accumulation_from_direction" << std::endl;
  auto accumulation_image = accumulation_from_direction(model.direction());
  std::cout << "done it" << std::endl;
  Z accumulation_vec(gfs);
  ZDGF accumulation_dgf(gfs,accumulation_vec);
  for (const auto& e : elements(gv))
    {
      auto x = e.geometry().center();
      int cellx = std::floor(x[0]/90.0);
      int celly = std::floor(x[1]/90.0);
      Dune::PDELab::Backend::native(accumulation_vec)[indexset.index(e)][0] = accumulation_image(cellx,celly);
    }

  // pixels upstream of outlet
  std::vector<std::pair<int,int>> outlet_pixels;
  {
    std::pair<int,int> pixel;
    pixel.first = std::floor(15075.0/90.0);
    pixel.second = std::floor(16875.0/90.0);
    int nprobes=5;
    outlet_pixels.push_back(pixel);
    for (int i=1; i<=nprobes; i++)
      {
        std::pair<int,int> pixel2(pixel);
        pixel2.second += i;
        outlet_pixels.push_back(pixel2);
        pixel2 = pixel;
        pixel2.second -= i;
        outlet_pixels.push_back(pixel2);
      }
  }
  std::cout << "call upstream_from_direction" << std::endl;
  auto upstream_image = upstream_from_direction(model.direction(),outlet_pixels);
  Z upstream_vec(gfs);
  ZDGF upstream_dgf(gfs,upstream_vec);
  for (const auto& e : elements(gv))
    {
      auto x = e.geometry().center();
      int cellx = std::floor(x[0]/90.0);
      int celly = std::floor(x[1]/90.0);
      Dune::PDELab::Backend::native(upstream_vec)[indexset.index(e)][0] = upstream_image(cellx,celly);
    }

  // curvature
  auto curvature_image = curvature(model.elevation());
  Z curvature_vec(gfs);
  ZDGF curvature_dgf(gfs,curvature_vec);
  for (const auto& e : elements(gv))
    {
      auto x = e.geometry().center();
      int cellx = std::floor(x[0]/90.0);
      int celly = std::floor(x[1]/90.0);
      Dune::PDELab::Backend::native(curvature_vec)[indexset.index(e)][0] = curvature_image(cellx,celly);
    }

  // gradient
  auto flatness_image = flatness(model.elevation());
  Z flatness_vec(gfs);
  ZDGF flatness_dgf(gfs,flatness_vec);
  for (const auto& e : elements(gv))
    {
      auto x = e.geometry().center();
      int cellx = std::floor(x[0]/90.0);
      int celly = std::floor(x[1]/90.0);
      Dune::PDELab::Backend::native(flatness_vec)[indexset.index(e)][0] = flatness_image(cellx,celly);
    }

  // snow cover
  Z snowcover_vec(gfs);
  ZDGF snowcover_dgf(gfs,snowcover_vec);
  for (const auto& e : elements(gv))
    {
      auto x = e.geometry().center();
      int cellx = std::floor(x[0]/90.0);
      int celly = std::floor(x[1]/90.0);
      Dune::PDELab::Backend::native(snowcover_vec)[indexset.index(e)][0] = model.snowcover(0)(cellx,celly);
    }

  // glacier cover
  Z glacier_vec(gfs);
  ZDGF glacier_dgf(gfs,glacier_vec);
  for (const auto& e : elements(gv))
    {
      auto x = e.geometry().center();
      int cellx = std::floor(x[0]/90.0);
      int celly = std::floor(x[1]/90.0);
      Dune::PDELab::Backend::native(glacier_vec)[indexset.index(e)][0] = model.glacier()(cellx,celly);
    }

  // temperature
  Z temperature_vec(gfs);
  ZDGF temperature_dgf(gfs,temperature_vec);
  for (const auto& e : elements(gv))
    {
      auto x = e.geometry().center();
      int cellx = std::floor(x[0]/90.0);
      int celly = std::floor(x[1]/90.0);
      Dune::PDELab::Backend::native(temperature_vec)[indexset.index(e)][0] = model.temperature(0)(cellx,celly);
    }

  // precipitation
  Z precipitation_vec(gfs);
  ZDGF precipitation_dgf(gfs,precipitation_vec);
  for (const auto& e : elements(gv))
    {
      auto x = e.geometry().center();
      int cellx = std::floor(x[0]/90.0);
      int celly = std::floor(x[1]/90.0);
      Dune::PDELab::Backend::native(precipitation_vec)[indexset.index(e)][0] = model.precipitation(0)(cellx,celly);
    }

  // lst
  Z lst_vec(gfs);
  ZDGF lst_dgf(gfs,lst_vec);
  for (const auto& e : elements(gv))
    {
      auto x = e.geometry().center();
      int cellx = std::floor(x[0]/90.0);
      int celly = std::floor(x[1]/90.0);
      Dune::PDELab::Backend::native(lst_vec)[indexset.index(e)][0] = model.land_surface_temperature(0)(cellx,celly);
    }

  // source term
  Z source_vec(gfs);
  ZDGF source_dgf(gfs,source_vec);
  Dune::FieldVector<double,2> x0(0.5);
  for (const auto& e : elements(gv))
    Dune::PDELab::Backend::native(source_vec)[indexset.index(e)][0] = model.f(e,x0);

  // coordinates of pixel
  Z long_vec(gfs);
  ZDGF long_dgf(gfs,long_vec);
  Z lat_vec(gfs);
  ZDGF lat_dgf(gfs,lat_vec);
  double ox = ptree.get("problem.ox",(double)77.3987); 
  double oy = ptree.get("problem.oy",(double)33.89998); 
  double hx = 90.0;
  double hy = 90.0;
  for (const auto& e : elements(gv))
    {
      auto x = e.geometry().center();
      Dune::PDELab::Backend::native(long_vec)[indexset.index(e)][0] = ox+(x[0]/hx/1200.0);
      Dune::PDELab::Backend::native(lat_vec)[indexset.index(e)][0] = oy+(x[1]/hy/1200.0);
    }
  
  // new writer
  std::cout << "make vtk writer" << std::endl;
  using Writer = Dune::VtkImageDataWriter<GV>;
  Dune::PvdWriter<Writer> pvdWriter(gv, Dune::Vtk::FormatTypes::COMPRESSED, Dune::Vtk::DataTypes::FLOAT32);
  pvdWriter.addCellData(std::make_shared<VTKF20>(z2newdgf0,"u_surface"));
  pvdWriter.addCellData(std::make_shared<VTKF21>(z2newdgf1,"u_groundwater"));
  pvdWriter.addCellData(std::make_shared<VTKF>(long_dgf,"longitude e"));
  pvdWriter.addCellData(std::make_shared<VTKF>(lat_dgf,"latitude n"));
  pvdWriter.addCellData(std::make_shared<VTKF>(bsurfacedgf,"bathymmetry_surface"));
  pvdWriter.addCellData(std::make_shared<VTKF>(bgroundwaterdgf,"bathymmetry_groundwater"));
  pvdWriter.addCellData(std::make_shared<VTKF>(hdgf_surface,"h_surface"));
  pvdWriter.addCellData(std::make_shared<VTKF>(hdgf_groundwater,"h_groundwater"));
  pvdWriter.addCellData(std::make_shared<VTKF>(catchment_dgf,"catchment"));
  pvdWriter.addCellData(std::make_shared<VTKF>(snowcover_dgf,"snow cover"));
  pvdWriter.addCellData(std::make_shared<VTKF>(glacier_dgf,"glacier"));
  pvdWriter.addCellData(std::make_shared<VTKF>(temperature_dgf,"temperature"));
  pvdWriter.addCellData(std::make_shared<VTKF>(lst_dgf,"land surface temperature"));
  pvdWriter.addCellData(std::make_shared<VTKF>(precipitation_dgf,"precipitation"));
  pvdWriter.addCellData(std::make_shared<VTKF>(source_dgf,"source"));
  pvdWriter.addCellData(std::make_shared<VTKF>(accumulation_dgf,"accumulation"));
  pvdWriter.addCellData(std::make_shared<VTKF>(upstream_dgf,"upstream"));
  pvdWriter.addCellData(std::make_shared<VTKF>(curvature_dgf,"curvature"));
  pvdWriter.addCellData(std::make_shared<VTKF>(flatness_dgf,"flatness"));
  pvdWriter.addCellData(std::make_shared<VeloVTKFSurface>(velodgfsurface,"velocity_surface"));
  pvdWriter.addCellData(std::make_shared<VeloVTKFGroundwater>(velodgfgroundwater,"velocity_groundwater"));
  std::string fullfilename = filename + ".vti";
  bool outputon = ptree.get<bool>("output.on");
  if (outputon)
    pvdWriter.writeTimestep(time,fullfilename);
  
  // time loop
  std::cout << "start time loop" << std::endl;
  RF T = ptree.get("problem.T",(RF)1.0);
  int simulation_days = ptree.get<int>("problem.simulation_days");
  T = simulation_days*86400.0;
  RF dt = ptree.get("method.dt",(RF)0.1);
  RF timestepmax =  ptree.get("method.dtmax",dt);
  int every = ptree.get("output.every",(int)1);
  int maxsteps = ptree.get("problem.maxsteps",(int)10000000);
  int step=0;
  int increased_step = 0;
  int decreased_step = 0;
  double eps=1e-6;
  auto factor = 1.0/std::sqrt(2.0);
  double accumulated_discharge_surface = 0.0;
  double accumulated_discharge_groundwater = 0.0;
  int next_day_boundary = (time/86400.0)+1;

  std::vector<std::pair<double,double>> simulated;
  simulated.push_back(std::pair<double,double>(time,0.0));
    
  while (time<T-1e-8 && step<maxsteps)
    {
      // adjust time step to hit day boundary
      // what is the next date boundary to hit?
      double time_to_hit = next_day_boundary*86400.0;

      if (time+dt>time_to_hit-1e-5)
        {
          dt = time_to_hit-time;
        }
      else if (time+2*dt>time_to_hit-1e-5)
        {
          dt = 0.5*(time_to_hit-time);
        }        
      
      // assemble constraints for new time step
      model.setTime(time+dt);

      // do time step
      try {
        z2new = z2;
        osm.apply(time,dt,z2,z2new);
      }
      catch (Dune::PDELab::LineSearchError) {
        dt *= factor;
        if (dt<eps) throw;
        decreased_step = step;
        if (gv.comm().rank()==0)
          std::cout << "STEP CONTROL: " << " decreased in step " << step+1 << " newdt: " << dt << std::endl;
        continue;
      }
      catch (Dune::PDELab::TerminateError) {
        dt *= factor;
        if (dt<eps) throw;
        decreased_step = step;
        if (gv.comm().rank()==0)
          std::cout << "STEP CONTROL: " << " decreased in step " << step+1 << " newdt: " << dt << std::endl;
        continue;
      }
      catch (Dune::PDELab::NewtonLinearSolverError) {
        dt *= factor;
        if (dt<eps) throw;
        decreased_step = step;
        if (gv.comm().rank()==0)
          std::cout << "STEP CONTROL: " << " decreased in step " << step+1 << " newdt: " << dt << std::endl;
        continue;
      }
      catch (Dune::ISTLError) {
        dt *= factor;
        if (dt<eps) throw;
        decreased_step = step;
        if (gv.comm().rank()==0)
          std::cout << "STEP CONTROL: " << " decreased in step " << step+1 << " newdt: " << dt << std::endl;
        continue;
      }

      // accept time step
      time+=dt;
      step++;
      z2 = z2new;
      if (std::abs(time-next_day_boundary*86400.0)<1e-4)
        next_day_boundary += 1; // we have finished one day

      // split coefficient vector
      Dune::PDELab::interpolate(z2newdgf0,gfs,z_surface);
      Dune::PDELab::interpolate(z2newdgf1,gfs,z_groundwater);

      // update velocity fields
      velodgfsurface.update();
      velodgfgroundwater.update();

      // compute height
      h_surface = z_surface;
      h_surface -= bathymmetry_surface;
      h_groundwater = z_groundwater;
      h_groundwater -= bathymmetry_groundwater;

      // update temperature and snow cover
      int timestep = day_from_time(time);
      timestep = std::min(timestep,637);

      for (const auto& e : elements(gv))
        {
          auto x = e.geometry().center();
          int cellx = std::floor(x[0]/90.0);
          int celly = std::floor(x[1]/90.0);
          auto index = indexset.index(e);
          Dune::PDELab::Backend::native(snowcover_vec)[index][0] = model.snowcover(timestep)(cellx,celly);
          Dune::PDELab::Backend::native(temperature_vec)[index][0] = model.temperature(timestep)(cellx,celly);
          Dune::PDELab::Backend::native(precipitation_vec)[index][0] = model.precipitation(timestep)(cellx,celly);
          Dune::PDELab::Backend::native(lst_vec)[index][0] = model.land_surface_temperature(timestep)(cellx,celly);
          Dune::PDELab::Backend::native(source_vec)[index][0] = model.f_surface(e,x0);
        }

      // probe for discharge output
      using VeloProbeSurface = Dune::PDELab::GridFunctionProbe<VeloDGFSurface>;
      using VeloProbeGroundwater = Dune::PDELab::GridFunctionProbe<VeloDGFGroundwater>;
      using Domain = typename VeloDGFSurface::Traits::DomainType;
      using VeloRange = typename VeloDGFSurface::Traits::RangeType;
      std::vector<VeloProbeSurface> probevec_surface;
      std::vector<VeloProbeGroundwater> probevec_groundwater;
      Domain point;
      point[0] = 15075.0; point[1] = 16875.0; // the discharge measurement point
      int nprobes=5; // in each direction up and down
      probevec_surface.push_back(VeloProbeSurface(velodgfsurface,point));
      probevec_groundwater.push_back(VeloProbeGroundwater(velodgfgroundwater,point));
      for (int i=1; i<=nprobes; i++)
        {
          Domain point2(point);
          point2[1] += i*90.0;
          probevec_surface.push_back(VeloProbeSurface(velodgfsurface,point2));
          probevec_groundwater.push_back(VeloProbeGroundwater(velodgfgroundwater,point2));
          point2 = point;
          point2[1] -= i*90.0;
          probevec_surface.push_back(VeloProbeSurface(velodgfsurface,point2));
          probevec_groundwater.push_back(VeloProbeGroundwater(velodgfgroundwater,point2));
        }
      RF discharge_surface = 0.0;
      RF discharge_groundwater = 0.0;
      for (auto& probe : probevec_surface)
        {
          VeloRange result(0.0);
          probe.eval_all(result);
          discharge_surface += result[0]*90.0;
        }
      for (auto& probe : probevec_groundwater)
        {
          VeloRange result(0.0);
          probe.eval_all(result);
          discharge_groundwater += result[0]*90.0;
        }
      accumulated_discharge_surface += dt*discharge_surface;
      accumulated_discharge_groundwater += dt*discharge_groundwater;
      if (time/86400.0>365.0 && time/86400.0<375.0) accumulated_discharge_surface = 0.0;
      if (time/86400.0>365.0 && time/86400.0<375.0) accumulated_discharge_groundwater = 0.0;

      // compute water stored in the catchment
      RF stored_surfacewater=0.0;
      RF stored_groundwater=0.0;
      for (const auto& e : elements(gv,Dune::Partitions::interior))
        {
          auto x = e.geometry().center();
          int cellx = std::floor(x[0]/90.0);
          int celly = std::floor(x[1]/90.0);
          auto index = indexset.index(e);
          if (Dune::PDELab::Backend::native(upstream_vec)[index][0]!=1) continue;
          auto h_surface_cell = Dune::PDELab::Backend::native(h_surface)[index][0];
          auto h_groundwater_cell = Dune::PDELab::Backend::native(h_groundwater)[index][0];
          Domain xlocal(0.5);
          auto storage_surface_cell = model.storage_surface(e,xlocal);
          auto storage_groundwater_cell = model.storage_groundwater(e,xlocal);
          stored_surfacewater += storage_surface_cell*h_surface_cell*90.0*90.0;
          stored_groundwater += storage_groundwater_cell*h_groundwater_cell*90.0*90.0;
        }
      stored_surfacewater = gv.comm().sum(stored_surfacewater);
      stored_groundwater = gv.comm().sum(stored_groundwater);
      
      
      // output discharge and accumulated fluids
      if (gv.comm().rank()==0)
        std::cout << "DISCHARGECOUPLED "
                  << time/86400.0
                  << " " << discharge_surface
                  << " " << discharge_groundwater
                  << " " << accumulated_discharge_surface
                  << " " << accumulated_discharge_groundwater
                  << " " << stored_surfacewater
                  << " " << stored_groundwater
                  << std::endl;
      simulated.push_back(std::pair<double,double>(time,discharge_surface));
      
      // output to VTK file
      if (outputon && step%every==0)
        {
          pvdWriter.writeTimestep(time,fullfilename);
          if (gv.comm().rank()==0)
            std::cout << "WRITING VTK OUTPUT " << step << " " << time << " " << dt << std::endl;
        }

      // increase timestep
      if (step-decreased_step>5 && step-increased_step>5)
        {
          double newdt = std::min(timestepmax,std::sqrt(2.0)*dt);
          if (newdt>dt)
            {
              increased_step = step;
              if (gv.comm().rank()==0)
                std::cout << "STEP CONTROL: " << " increased in step " << step << " newdt: " << newdt << std::endl;
              dt = newdt;
            }
        }
    }

  // write state to file
  std::string storefilename=ptree.get<std::string>("output.storefilename");
  auto s = storefilename + "_" + std::to_string(gv.comm().rank()) + ".mm";
  storeMatrixMarket(Dune::PDELab::Backend::native(z2new),s);

  // prepare return value
  DataAnalysis analysis;
  std::map<int,double> averages = analysis.daily_averages(simulated);
  return averages;
}

std::string double_to_string (double x)
{
	char buffer[256];
	std::sprintf(buffer,"%.15e",x);
	return std::string(buffer);
}

// class that does a full simulation run of the coupled simulation
// for a specific set of parameters
template<typename Real, typename Point>
class Simulation
{
  std::vector<std::string> param_name;

public:

  // need to know the parameter names that correspond to the parameter values
  Simulation (std::vector<std::string> names)
    : param_name(names)
  {
  }

  // run with set of parameters given a tag
  Real run (Point p, std::string tag) const
  {
    // check correct number of parameters
    if (p.size()!=param_name.size()) return -1.0;
   
    // load parameter tree with baseline parameters
    Dune::ParameterTree ptree;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("stok.ini",ptree);

    // get basename
    std::string basename=ptree.get<std::string>("output.filename");

    // map parameters into parameter tree
    for (int i=0; i<p.size(); ++i)
      ptree[param_name[i]] = double_to_string(p[i]);

    // do a simulation of the first year to determine starting value
    ptree["problem.tstart"] = "0.0";
    ptree["problem.simulation_days"] = "365";
    ptree["output.loadinitial"] = "false";
    ptree["output.storefilename"] = basename + "_" + tag + "_" + "warmup";
    auto a_warmup = driverCoupledFVNewNewton(ptree);

    // then do the full run
    ptree["problem.tstart"] = "0.0";
    ptree["problem.simulation_days"] = "638";
    ptree["output.loadinitial"] = "true";
    ptree["output.loadfilename"] = basename + "_" + tag + "_" + "warmup";
    ptree["output.storefilename"] = basename + "_" + tag + "_" + "truerun";
    auto a_truerun = driverCoupledFVNewNewton(ptree);

    // analysis
    std::string path_to_data = ptree.get<std::string>("problem.path_to_data");
    std::string data_file_name = ptree.get<std::string>("problem.data_file_name");
    DataAnalysis analysis;
    auto m = analysis.load_measurement_file(path_to_data + data_file_name);
    auto objective = analysis.l2_distance(m,a_truerun,365);

    // write averages
    auto averages_file_name = basename + "_" + tag + "_" + "averages.dat";
    std::ofstream myfile0;
    myfile0.open(averages_file_name);
    for (auto it=a_truerun.begin(); it!=a_truerun.end(); ++it)
      myfile0 << it->first << " " << it->second << std::endl;
    myfile0.close();

    // write objective value to file
    auto objective_file_name = basename + "_" + tag + "_" + "objective.dat";
    std::ofstream myfile1;
    myfile1.open(objective_file_name);
    myfile1 << objective
      << " p0=" << p[0]
      << " p1=" << p[1]
      << " p2=" << p[2]
      << std::endl;
    myfile1.close();

    // return result
    return objective;
  }
  
};

// problem class for nonlinopt
template<typename Real, typename Point>
class Problem
  : public Dune::NonlinOpt::FiniteDifferenceProblemBase<Real,Point>
{
  std::string tag;
  Point lows,highs;
  int rank;
  Simulation<Real,Point> simulation;
  mutable int counter;

public:
  Problem (std::string tag_, std::vector<std::string> names, Point lows_, Point highs_, int rank_)
    : Dune::NonlinOpt::FiniteDifferenceProblemBase<Real,Point>({1e-2,1e-2,1e-2}), tag(tag_), lows(lows_), highs(highs_), rank(rank_), 
        simulation(names), counter(1)
  {}
  
  Real value(const Point& point, bool subsequent = false) const override
  {
    // check if point is in range
    // if not, return infinity
    for (int i=0; i<point.size(); ++i)
      if (point[i]<lows[i] || point[i]>highs[i])
        return std::numeric_limits<double>::infinity();

    // assume 2nd parameter is log conductivity
    Point pt(point);
    pt[2] = std::exp(pt[2]);
        
    // now we can do a simulation
    try {
      if (rank==0)
        std::cout << "RUNNING SIMULATION: p0=" << pt[0] << " p1=" << pt[1] << " p2=" << pt[2] << std::endl;
      auto J = simulation.run(pt,tag+std::to_string(counter));
      counter++;
      return J;
    }
    // which still might fail
    catch (...)
      {
        return std::numeric_limits<double>::infinity();
      }
  }

  Point zero() const override
  {
    return Point(dim());
  }
  
  std::size_t dim() const override
  {
    return 3;
  }
  
  void hook(std::size_t iteration, const Point& point,
            Real value, const Point& gradient, bool extrapolation = false) const override
  {
    if (rank==0)
      std::cout << "NONLINOPT iter=" << iteration
                << " value=" << value
                << " p0=" << point[0]
                << " p1=" << point[1]
                << " p2=" << point[2]
                << std::endl;
  }
};

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
    ptreeparser.readINITree("stok.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);
    std::string mode=ptree.get<std::string>("problem.mode");

    // compute objective
    if (mode=="analysis")
      {
        std::string path_to_data = ptree.get<std::string>("problem.path_to_data");
        std::string data_file_name = ptree.get<std::string>("problem.data_file_name");
        DataAnalysis analysis;
        auto m = analysis.load_measurement_file(path_to_data + data_file_name);
        // for (auto x : m)
        //   {
        //     std::cout << x.first-1 << " " << x.second << std::endl;
        //     std::cout << x.first << " " << x.second << std::endl;
        //   }
        std::string simulation_file_name = ptree.get<std::string>("analysis.simulation_file_name");
        auto s = analysis.load_data_file(simulation_file_name);
        auto a = analysis.daily_averages(s);
        for (auto x : a)
          {
            std::cout << x.first-1 << " " << x.second << std::endl;
            std::cout << x.first << " " << x.second << std::endl;
          }
        std::cout << "OBJECTIVE from averages " << analysis.l2_distance(m,a) << std::endl;
        return 0;
      }
    
    // write parameter tree to output
    if (helper.rank()==0)
      {
        std::cout << "******************** parameter tree ********************" << std::endl;
        ptree.report();
        std::cout << "********************************************************" << std::endl;
      }

    // read mode and call simulation
    if (mode=="scalar")
      {
        driverFVNewNewton(ptree);
        return 0;
      }
    if (mode=="coupled")
      {
        using Point = std::vector<double>;
        Point p;
        p.push_back(ptree.get<double>("problem.k_surface"));
        p.push_back(ptree.get<double>("problem.d0"));
        p.push_back(ptree.get<double>("problem.L_i"));
        std::vector<std::string> param_names = {"problem.k_surface","problem.d0","problem.L_i"};
        Simulation<double,Point> simulation(param_names);
        std::string tag=ptree.get<std::string>("output.tag");
        auto J = simulation.run(p,tag);
        std::cout << "J=" << J << std::endl;
        return 0;
      }
    
    if (mode=="nonlinopt")
      {
        // type for parameter value
        using Real = double;
        using Point = Dune::NonlinOpt::VectorClass<Real>;

        // initial value
        Point point = {ptree.get<double>("problem.k_surface"),
                       ptree.get<double>("problem.d0"),
                       std::log(ptree.get<double>("problem.L_i"))
        };

        // parameter ranges
        Point lows = {0.1,4.0,-30.0};
        Point highs = {10.0,40.0,-7.0};

	// parameter names
        std::vector<std::string> param_names = {"problem.k_surface","problem.d0","problem.L_i"};

        // problem class instance
        std::string tag=ptree.get<std::string>("output.tag");
        Problem<Real,Point> problem(tag,param_names,lows,highs,helper.rank());

        // make solver
        using Solver = Dune::NonlinOpt::UnconstrainedOptimization<Real,Point>;
        Solver solver(ptree);
        solver.report();

        // call solver
        solver.apply(problem,point);
        
        return 0;
      }

    // std::vector<std::vector<double>> paramsA = {{0.5,1.0,1.3},{10.0,12.6,16.0},{4e-5,1e-5,0.25e-5}};
    std::vector<std::vector<double>> params = {{1.6,1.8,2.0,2.2,2.4},{17.5,20.0,22.5},{1.6e-6,1.8e-6,2e-6,2.2e-6,2.4e-6}};
    std::vector<std::string> param_names = {"problem.k_surface","problem.d0","problem.L_i"};
    if (mode=="parameter_grid")
      {

        int p0 = ptree.get<int>("param.p0"); if (p0<0 || p0>=params[0].size()) exit(1);
        int p1 = ptree.get<int>("param.p1"); if (p1<0 || p1>=params[1].size()) exit(1);
        int p2 = ptree.get<int>("param.p2"); if (p2<0 || p2>=params[2].size()) exit(1);
        using Point = std::vector<double>;
        Point p;
        p.push_back(params[0][p0]);
        p.push_back(params[1][p1]);
        p.push_back(params[2][p2]);
        std::string tag = std::to_string(p0) + std::to_string(p1) + std::to_string(p2);
        Simulation<double,Point> simulation(param_names);
        std::cout << "Running simulation tag=" << tag
                  << " " << param_names[0] << "=" << p[0]
                  << " " << param_names[1] << "=" << p[1]
                  << " " << param_names[2] << "=" << p[2]
                  << std::endl;
        auto J = simulation.run(p,tag);
        std::cout << "J=" << J << std::endl;
        return 0;
      }
    if (mode=="parameter_script")
      {
        std::string basename=ptree.get<std::string>("output.filename");
        for (int i0=0; i0<params[0].size(); i0++)
          for (int i1=0; i1<params[1].size(); i1++)
            for (int i2=0; i2<params[2].size(); i2++)
              {
                std::string tag = std::to_string(i0) + std::to_string(i1) + std::to_string(i2);
                std::cout << "nohup ./stok"
                          << " -param.p0 " << i0 
                          << " -param.p1 " << i1
                          << " -param.p2 " << i2
                          << " > "
                          << basename + "_" + tag + ".log &"
                          << std::endl;
              }
      }
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
