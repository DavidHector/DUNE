// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_HYDRO_LOCALOPERATOR_NONLINEARDIFFUSIONCOUPLEDMODEL_HH
#define DUNE_HYDRO_LOCALOPERATOR_NONLINEARDIFFUSIONCOUPLEDMODEL_HH

#include<dune/common/parametertreeparser.hh>

/** a class that can be passed as parameter to NonlinearDiffusionFVCoupled
 *
 * This example solves the same problem in groundwater and surface
 * it also offers the scalar functions so that it can be passed to NonlinearDiffusionFV
 * Thus it can be used to test the scalar and coupled problems
 */
template<typename Number>
class NonlinearDiffusionCoupledModel
{
  Number eps,A,B;

public:
  //! Constructor without arg sets nonlinear term to zero
  NonlinearDiffusionCoupledModel () : eps(1e-3)
  {
    // regularization parameter 
    A = 2.5/(2.0*std::sqrt(eps));
    B = 0.5/(2.0*pow(eps,2.5));
    std::cout << "REGULARIZATION gradient eps=" << eps << " A=" << A << " B=" << B << std::endl;
  }

  // We provide functions for both equations as before plus the coupling term
  
  //*****************************************************************************
  // surface equation
  //*****************************************************************************

  //! storage term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number storage_surface (const E& e, const X& x) const
  {
    return 1.0;
  }

  //! bottom position in global coordinates (function b above)
  /**
     \param e      : codim 0 entity
     \param xlocal : local coordinate
   */
  template<typename E, typename X>
  Number bathymmetry_surface (const E& e, const X& xlocal) const
  {
    auto x = e.geometry().center();
    return 5.0 + x[0]*0.01 +  (x[1]-25.0)*(x[1]-25.0)/625.0*50.0 ;
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
    // upwind evaluation of height
    Number u_upwind = (vn>=0) ? u_inside : u_outside;
    
    // evaluation of (u-b)
    Number bmax=std::max(b_inside,b_outside);
    Number h_upwind = std::max(u_upwind-bmax,0.0);
    
    // regularization of gradient and evaluation
    Number gradterm =  (std::abs(vn)>=eps) ? 1.0/std::sqrt(std::abs(vn)) : (A-B*vn*vn);

    return h_upwind*gradterm;
  }

  //! permeability value at position in a cell
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number k_surface (const E& e, const X& x) const
  {
    return 1.0;
  }
  
  //! source/sink term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number f_surface (const E& e, const X& x) const
  {
    auto center = e.geometry().center();
    return 0.0;
  }
  
  //! boundary condition type function (true = Dirichlet)
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  bool b_surface (const I& i, const X& x) const
  {
    return true;
  }

  //! Dirichlet extension
  /** Provides Dirichlet boundary values as well as initial value!
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number g_surface (const E& e, const X& x) const
  {
    auto center = e.geometry().center();
    auto z(center);
    z[0] = 175.0; z[1] = 25.0;
    z -= center;
    auto distance = z.two_norm();
    return bathymmetry_surface(e,x) + std::max(0.25*(1.0-distance/15.0),0.0);
  }

  //! Neumann boundary condition
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  Number j_surface (const I& i, const X& x) const
  {
    return 0.0;
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
    return storage_surface(e,x);
  }

  //! bottom position in global coordinates (function b above)
  /**
     \param e      : codim 0 entity
     \param xlocal : local coordinate
   */
  template<typename E, typename X>
  Number bathymmetry_groundwater (const E& e, const X& xlocal) const
  {
    return bathymmetry_surface(e,xlocal);
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
    return phi_surface(b_inside,b_outside,u_inside,u_outside,vn);
  }

  //! permeability value at position in a cell
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number k_groundwater (const E& e, const X& x) const
  {
    return k_surface(e,x);
  }
  
  //! source/sink term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number f_groundwater (const E& e, const X& x) const
  {
    return f_surface(e,x);
  }
  
  //! boundary condition type function (true = Dirichlet)
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  bool b_groundwater (const I& i, const X& x) const
  {
    return b_surface(i,x);
  }

  //! Dirichlet extension
  /** Provides Dirichlet boundary values as well as initial value!
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number g_groundwater (const E& e, const X& x) const
  {
    return g_surface(e,x);
  }

  //! Neumann boundary condition
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  Number j_groundwater (const I& i, const X& x) const
  {
    return j_surface(i,x);
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
    return storage_surface(e,x);
  }

  //! bottom position in global coordinates (function b above)
  /**
     \param e      : codim 0 entity
     \param xlocal : local coordinate
   */
  template<typename E, typename X>
  Number bathymmetry (const E& e, const X& xlocal) const
  {
    return bathymmetry_surface(e,xlocal);
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
    return phi_surface(b_inside,b_outside,u_inside,u_outside,vn);
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
    return k_surface(e,x);
  }
  
  //! source/sink term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number f (const E& e, const X& x) const
  {
    return f_surface(e,x);
  }
  
  //! boundary condition type function (true = Dirichlet)
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  bool b (const I& i, const X& x) const
  {
    return b_surface(i,x);
  }

  //! Dirichlet extension
  /** Provides Dirichlet boundary values as well as initial value!
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number g (const E& e, const X& x) const
  {
    return g_surface(e,x);
  }

  //! Neumann boundary condition
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  Number j (const I& i, const X& x) const
  {
    return j_surface(i,x);
  }

  //*****************************************************************************
  // and the coupling
  //*****************************************************************************

  // Exchange term
  template<typename E, typename X>
  Number q (const E& e, const X& x, Number u_surface, Number u_groundwater, Number bath_surf) const
  {
    return 0.0;
  }

  void setTime (double t)
  {}
};

#endif
