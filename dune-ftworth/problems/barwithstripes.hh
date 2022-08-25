#ifndef Udune_ftworth_barwithstripes_HH
#define Udune_ftworth_barwithstripes_HH

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
class BarWithStripes
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  BarWithStripes(Dune::ParameterTree& ptree_) {}

  //! tensor diffusion constant per cell? return false if you want more than one evaluation of A per cell.
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::PermTensorType I;
    typename Traits::DomainType xglobal = e.geometry().center();
    // initialize isotropic tensor from field
    //double k1 = 1;
    //double k2 = contrast;
    //yfg.evaluateK(xglobal,k1);
    typename Traits::PermTensorType K(0.0);
    I[0][0] = 1.0;
    I[0][1] = 0.0;
    I[1][0] = 0.0;
    I[1][1] = 1.0;

    // double l0 = 8.0; // extension of domain in x0-direction
    // double l1 = 1.0; // extension of domain in x1-direction
    double contrast = 10000.0;

    // put some high permeability regions into the overlap

    // horizontal stripes
    for (int i = 0; i<10; i++) {
      if ( xglobal[1] < i*0.1+0.05 && xglobal[1]> i*0.1) {
          I[0][0] = contrast;
          I[1][1] = contrast;
          return I;         
      }
    }

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
    typename Traits::DomainType xglobal = e.geometry().center();
    if (xglobal[1] > 0.6)
      return 1000.0;
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
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    typename Traits::DomainType xglobal = is.geometry().global(x);
    if (xglobal[0]<1e-6) return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    if (xglobal[0]>8.0-1e-6) return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(xlocal);
    return xglobal[0];
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
