// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_HYDRO_LOCALOPERATOR_NONLINEARDIFFUSIONMODEL_HH
#define DUNE_HYDRO_LOCALOPERATOR_NONLINEARDIFFUSIONMODEL_HH

/** a class that can be passed as parameter to NonlinearDiffusionFV
 *
 * Provides the parameters for
 * \f{align*}{
 * \partial_t(s(x) u) - \nabla\cdot(k(x)\varphi(x,u,\nabla u)\nabla u) &= f(x,t) &&\text{in $\Omega\times\Sigma$}, \label{eq:flow}\ \
 * u(x,t_0) &= u_0(x) &&\text{at $t=t_0$}, \label{eq:ic} \\
 * u(x,t) &= g(x,t) &&\text{for $x\in\Gamma_D(t)$}, \label{eq:dirichlet}\\
 * - (k(x)\varphi(x,u,\nabla u)\nabla u)\cdot\nu &= j(x,t) &&\text{for $x\in\Gamma_N(t)$}, \label{eq:neumann}
 * \f}
 *
 */

// define the problem
template<typename Number>
class Model
{
public:

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
  /**
     \param e      : codim 0 entity
     \param xlocal : local coordinate
   */
  template<typename E, typename X>
  Number bathymmetry (const E& e, const X& xlocal) const
  {
    return 0.0;
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
    
    return h_upwind;
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
    
    return 1.0;
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

    return 1.0;
  }

  template<typename U, typename VN>
  Number dphi_dvn (U b_inside, U b_outside, U u_inside, U u_outside, VN vn) const
  {
    return 0.0;
  }

  //! permeability value at position in a cell
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number k (const E& e, const X& x) const
  {
    return 1e-7;
  }
  
  //! source/sink term
  /**
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number f (const E& e, const X& x) const
  {
    return 0.0;
  }
  
  //! boundary condition type function (true = Dirichlet)
  /**
     \param i intersection (with boundary)
     \param x local coordinate on intersection
   */
  template<typename I, typename X>
  bool b (const I& i, const X& x) const
  {
    return false;
  }

  //! Dirichlet extension
  /** Provides Dirichlet boundary values as well as initial value!
     \param e codim 0 entity
     \param x local coordinate
   */
  template<typename E, typename X>
  Number g (const E& e, const X& x) const
  {
    return 10.0;
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

  //! set time for subsequent evaluation
  /**
     \param t value of time for subsequent evaluations
   */
  void setTime (double t)
  {}
};

#endif
