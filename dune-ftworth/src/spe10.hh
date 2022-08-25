#ifndef Udune_ftworth_spe10problem_HH
#define Udune_ftworth_spe10problem_HH

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

/** Parameter class for the stationary convection-diffusion equation of the following form:
 *
 * \f{align*}{
 *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \ \
 *                                              u &=& g \mbox{ on } \partial\Omega_D (Dirichlet)\ \
 *                (b(x,u) - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N (Flux)\ \
 *                        -(A(x)\nabla u) \cdot n &=& o \mbox{ on } \partial\Omega_O (Outflow)
 * \f}
 * Note:
 *  - This formulation is valid for velocity fields which are non-divergence free.
 *  - Outflow boundary conditions should only be set on the outflow boundary
 *
 * The template parameters are:
 *  - GV a model of a GridView
 *  - RF numeric type to represent results
 */

template<typename GV, typename RF>
class SPE10Problem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

  RF GRID_LX,GRID_LY,GRID_LZ;
  RF LX,LY,LZ;
  RF hx,hy,hz;
  RF gradx,grady;
  RF scaling;

  const int nx=60, ny=220, nz=85; // hx = 1200/60=20, hy=2200/220=10, hz=170/85=2
  const int N=nx*ny*nz;
  std::vector<float> kx,ky,kz;
  RF kxmean, kymean, kzmean;
  float kxmin, kymin, kzmin;
  float kxmax, kymax, kzmax;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  SPE10Problem (Dune::ParameterTree& ptree_) : ptree(ptree_), kx(N), ky(N), kz(N)
  {
    // works only in 3d
    if (Traits::dimDomain!=3) DUNE_THROW(Dune::Exception,"spe10: dim must be 3");

    // read structured grid parameters
    GRID_LX = ptree.get<double>("grid.LX");
    GRID_LY = ptree.get<double>("grid.LY");
    GRID_LZ = ptree.get<double>("grid.LZ");
    LX = ptree.get<double>("spe10.LX");
    LY = ptree.get<double>("spe10.LY");
    LZ = ptree.get<double>("spe10.LZ");
    hx = LX/nx;
    hy = LY/ny;
    hz = LZ/nz;

    std::string spe10filename = ptree.get<std::string>("spe10.filename");
    scaling = ptree.get<double>("spe10.scaling");

    FILE *ptr_file;
    ptr_file = fopen(spe10filename.c_str(),"r");
    if (!ptr_file) DUNE_THROW(Dune::Exception,"spe10: cannot open data file");

    std::cout << "start reading data file" << std::endl;
 
    // read kx
    kxmin=1e30; kxmax=-1e30;
    kxmean = 0.0;
    for (int k=0; k<nz; k++)
      for (int j=0; j<ny; j++)
	for (int i=0; i<nx; i++)
	  {
	    float perm;
	    if (fscanf(ptr_file,"%g",&perm)!=EOF)
	      {
		kx[(nz-1-k)*nx*ny+j*nx+i] = perm;
	      }
	    else
	      {
		std::cout << "kx: unexpected EOF i=" << i << " j=" << j << " k=" << k << std::endl;
		DUNE_THROW(Dune::Exception,"spe10: unexpected end of file");
	      }
	    kxmin=std::min(kxmin,perm);
	    kxmax=std::max(kxmax,perm);
	    kxmean += perm;
	  }
    std::cout << "kxmin=" << kxmin << " kxmax=" << kxmax << std::endl;

    // read ky
    kymin=1e30; kymax=-1e30;
    kymean = 0.0;
    for (int k=0; k<nz; k++)
      for (int j=0; j<ny; j++)
    	for (int i=0; i<nx; i++)
    	  {
    	    float perm;
    	    if (fscanf(ptr_file,"%g",&perm)!=EOF)
	      ky[(nz-1-k)*nx*ny+j*nx+i] = perm;
    	    else
    	      {
    		std::cout << "ky: unexpected EOF i=" << i << " j=" << j << " k=" << k << std::endl;
    		DUNE_THROW(Dune::Exception,"spe10: unexpected end of file");
    	      }
    	    kymin=std::min(kymin,perm);
    	    kymax=std::max(kymax,perm);
	    kymean += perm;
    	  }
    std::cout << "kymin=" << kymin << " kymax=" << kymax << std::endl;

    // read kz
    kzmin=1e30; kzmax=-1e30;
    kzmean = 0.0;
    for (int k=0; k<nz; k++)
      for (int j=0; j<ny; j++)
    	for (int i=0; i<nx; i++)
    	  {
    	    float perm;
    	    if (fscanf(ptr_file,"%g",&perm)!=EOF)
	      kz[(nz-1-k)*nx*ny+j*nx+i] = perm;
    	    else
    	      {
    		std::cout << "kz: unexpected EOF i=" << i << " j=" << j << " k=" << k << std::endl;
    		DUNE_THROW(Dune::Exception,"spe10: unexpected end of file");
    	      }
    	    kzmin=std::min(kzmin,perm);
    	    kzmax=std::max(kzmax,perm);
	    kzmean += perm;
    	  }
    std::cout << "kzmin=" << kzmin << " kzmax=" << kzmax << std::endl;

    kxmean /= N;
    kymean /= N;
    kzmean /= N;

    std::cout << "kxmean=" << kxmean << std::endl;
    std::cout << "kymean=" << kymean << std::endl;
    std::cout << "kzmean=" << kzmean << std::endl;
 
    fclose(ptr_file);
    std::cout << "read data file successfully" << std::endl; 

    // read regional gradient parameter
    gradx = ptree.get("gradient.gradx",(double)-1.0);
    grady = ptree.get("gradient.grady",(double) 0.0);
  }

  //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }


  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType center = e.geometry().center();
    int i = floor(center[0]/hx);
    int j = floor(center[1]/hy);
    int k = floor(center[2]/hz);
    if (i<0 || i>=nx || j<0 || j>=ny || k<0 || k>=nz)
      {
	std::cout << "access out of range: " << center << " i="  << i << " j=" << j << " k=" << k << std::endl;
	DUNE_THROW(Dune::Exception,"i,j,k out of range");
      }

    typename Traits::PermTensorType I(0.0);
    I[0][0] = (1.0-scaling)*kxmean + scaling*kx[k*nx*ny+j*nx+i];
    I[1][1] = (1.0-scaling)*kymean + scaling*ky[k*nx*ny+j*nx+i];
    I[2][2] = (1.0-scaling)*kzmean + scaling*kz[k*nx*ny+j*nx+i];

    // std::cout << center[0]/hx << " " << center[1]/hy << " " << center[2]/hz 
    // 	      << " i="  << i << " j=" << j << " k=" << k << " KX=" << I[0][0] << " KY=" << I[1][1] << " KZ=" << I[2][2] << std::endl;

    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType v(0.0);
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    return 0.0;
  }

  //! boundary condition type function
  /* return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet for Dirichlet boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann for flux boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow for outflow boundary conditions
   */
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    if (xglobal[1]<1e-6) return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    if (xglobal[1]>GRID_LY-1e-6) return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(xlocal);
    return gradx*xglobal[0] + grady*xglobal[1];
  }

  //! flux boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

private:
  Dune::ParameterTree& ptree;
};


#endif // Udune_ftworth_HH
