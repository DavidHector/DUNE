#ifndef Udune_ftworth_islandsproblem_HH
#define Udune_ftworth_islandsproblem_HH

#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>

/** Parameter class for the stationary convection-diffusion equation of the following form:
 *
 * \f{align*}{
 *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \ \
 *                                              u &=& g \mbox{ on } \partial\Omega_D (Dirichlet)\ \
 *                (b(x) u - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N (Flux)\ \
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
class IslandsProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

  double scaling;
  
public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  IslandsProblem (double scaling_) : scaling(scaling_)
  {}
  
  //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType xglobal = e.geometry().center();

    int ix = std::floor(15.0*xglobal[0]);
    int iy = std::floor(15.0*xglobal[1]);
    auto x = xglobal[0];
    auto y = xglobal[1];
    
    double kappa = 1.0;

    if (x>0.3 && x<0.9 && y>0.6-(x-0.3)/6 && y<0.8-(x-0.3)/6)
      kappa =  pow(10,5.0) * (x+y) * 10.0;

    if (x>0.1 && x<0.5 && y>0.1+x && y<0.25+x)
      kappa =  pow(10,5.0) * (1.0+7.0*y) ;

    if (x>0.5 && x<0.9 && y>0.15-(x-0.5)*0.25 && y<0.35-(x-0.5)*0.25)
      kappa =  pow(10,5.0)*2.5;

    if (ix%2==0 && iy%2==0)
      kappa = pow(10,5.0) * (1.0+ix+iy);

    kappa = kappa*scaling + (1.0-scaling);
    
    typename Traits::PermTensorType I;
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? kappa : 0;
    
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
    if (xglobal[0]<1e-6) return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    if (xglobal[0]>1.0-1e-6) return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(xlocal);
    return 1.0-xglobal[0];
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
};

#endif // Udune_ftworth_HH
