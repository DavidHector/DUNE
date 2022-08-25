// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_HYDRO_LOCALOPERATOR_NONLINEARDIFFUSIONFV_HH
#define DUNE_HYDRO_LOCALOPERATOR_NONLINEARDIFFUSIONFV_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/type.hh>

#include<dune/geometry/referenceelements.hh>
#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/variablefactories.hh>

/** a local operator for spatial part of the nonlinear diffusion equation with cell-centered FV
 *
 * \f{align*}{
 *   \nabla\cdot k(x) phi(x,u(x),\nabla u(x) \nabla u(x)) &=& f(x) x\in\Omega,  \\
 *                                u(x) &=& g(x) x\in\partial\Omega_D \\
 *  F(x,u(x),\nabla u(x)) \cdot \nu(x) &=& j(x) x\in\partial\Omega_N \\
 *
 * F(x,u,\nabla u(x)) = - k(x) (u-b)^\alpha |\nabla u(x)|^{\gamma-1} \nabla u(x)
 * \f}
 *
 */
template<typename Model, typename GFS>
class NonlinearDiffusionFV :
  // public Dune::PDELab::NumericalJacobianSkeleton<NonlinearDiffusionFV<Model,GFS> >,
  // public Dune::PDELab::NumericalJacobianBoundary<NonlinearDiffusionFV<Model,GFS> >,
  public Dune::PDELab::NumericalJacobianApplySkeleton<NonlinearDiffusionFV<Model,GFS> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<NonlinearDiffusionFV<Model,GFS> >,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::FullSkeletonPattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
  Model& model;        // parameter functions

  // evaluation of grid function from previous time step
  typedef typename GFS::Traits::GridView GV;
  typedef typename GV::IndexSet IndexSet;
  using Z = Dune::PDELab::Backend::Vector<GFS,double>;
  std::shared_ptr<const GFS> pgfs;
  const GV& gv;
  const IndexSet& indexset;
  Z bathymmetryvector;

public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  enum { doPatternSkeleton = true };

  // residual assembly flags
  enum { doLambdaVolume = true };
  enum { doLambdaBoundary = true };
  enum { doAlphaSkeleton  = true };
  enum { doAlphaBoundary  = true };

  //! constructor stores a reference to the parameter object
  NonlinearDiffusionFV (Model& model_, const GFS& gfs_)
    : model(model_), pgfs(stackobject_to_shared_ptr(gfs_)), gv(pgfs->gridView()), 
      indexset(gv.indexSet()), bathymmetryvector(*pgfs)
  {
    bathymmetryvector=0.0;
    auto bathymmetrylambda = [&](const auto& e, const auto& x)
      {return model.bathymmetry(e,x);};
    auto bathymmetrygf = Dune::PDELab::
      makeInstationaryGridFunctionFromCallable(gv,bathymmetrylambda,model);
    Dune::PDELab::interpolate(bathymmetrygf,*pgfs,bathymmetryvector);
  }

  // set the current time in parameter class
  void setTime (double t)
  {
    model.setTime(t);
  }

  //! residual contribution of volume source term
  template<typename EG, typename LFSV, typename R>
  void lambda_volume (const EG& eg, const LFSV& lfsv,
                      R& r) const
  {
    // center of reference element
    auto cellgeo = eg.geometry();
    auto cellcenterlocal =
      referenceElement(cellgeo).position(0,0);

    // accumulate residual
    auto f = model.f(eg.entity(),cellcenterlocal);
    r.accumulate(lfsv,0,-f*cellgeo.volume());
  }

  //! residual contribution of boundary integral (Neumann boundary condition)
  template<typename IG, typename LFSV, typename R>
  void lambda_boundary (const IG& ig, const LFSV& lfsv_i,
                        R& r_i) const
  {
    // face volume for integration
    auto facegeo = ig.geometry();
    auto facecenterlocal =
      referenceElement(facegeo).position(0,0);

    // evaluate boundary condition and quit on Dirichlet
    bool isdirichlet = model.b(ig.intersection(),facecenterlocal);

    if (isdirichlet) return;

    // contribution to residual from Neumann boundary
    auto j = model.j(ig.intersection(),facecenterlocal);
    r_i.accumulate(lfsv_i,0,j*facegeo.volume());
  }

  //! residual contribution from skeleton terms
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_skeleton (const IG& ig,
         const LFSU& lfsu_i, const X& x_i, const LFSV& lfsv_i,
         const LFSU& lfsu_o, const X& x_o, const LFSV& lfsv_o,
         R& r_i, R& r_o) const
  {
    // inside and outside cells
    auto cell_inside = ig.inside();
    auto cell_outside = ig.outside();

    // inside and outside cell geometries
    auto insidegeo = cell_inside.geometry();
    auto outsidegeo = cell_outside.geometry();

    // cell centers in local coordinates
    auto cell_center_local = referenceElement(insidegeo).position(0,0);
    
    // cell centers in global coordinates
    auto inside_global = insidegeo.center();
    auto outside_global = outsidegeo.center();

    // distance between the two cell centers
    inside_global -= outside_global;
    auto distance = inside_global.two_norm();

    // face volume for integration
    auto facegeo = ig.geometry();
    auto face_volume = facegeo.volume();
    
    // current solution
    auto u_inside = x_i(lfsu_i,0);
    auto u_outside = x_o(lfsu_o,0);

    // bathymmetry
    auto b_inside = Dune::PDELab::Backend::native(bathymmetryvector)[indexset.index(cell_inside)][0];
    auto b_outside = Dune::PDELab::Backend::native(bathymmetryvector)[indexset.index(cell_outside)][0];

    // gradient computed from current solution
    auto vn = (u_inside-u_outside)/distance;

    // permeability
    auto k_inside = model.k(cell_inside,cell_center_local);
    auto k_outside = model.k(cell_outside,cell_center_local);
    auto k = 2.0*(k_inside+1e-15)*(k_outside+1e-15)/(k_inside+k_outside+2e-15);
    
    // nonlinearity
    auto phi = model.phi(b_inside,b_outside,u_inside,u_outside,vn);
    
    // contribution to residual on inside and outside elements
    auto flux = k*phi*vn*face_volume;
    r_i.accumulate(lfsv_i,0, flux);
    r_o.accumulate(lfsv_o,0,-flux);
  }

  //! Jacobian contribution from skeleton terms
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename M>
  void jacobian_skeleton (const IG& ig,
         const LFSU& lfsu_i, const X& x_i, const LFSV& lfsv_i,
         const LFSU& lfsu_o, const X& x_o, const LFSV& lfsv_o,
         M& mat_ii, M& mat_io,
         M& mat_oi, M& mat_oo) const
  {
    // inside and outside cells
    auto cell_inside = ig.inside();
    auto cell_outside = ig.outside();

    // inside and outside cell geometries
    auto insidegeo = cell_inside.geometry();
    auto outsidegeo = cell_outside.geometry();

    // cell centers in local coordinates
    auto cell_center_local = referenceElement(insidegeo).position(0,0);
    
    // cell centers in global coordinates
    auto inside_global = insidegeo.center();
    auto outside_global = outsidegeo.center();

    // distance between the two cell centers
    inside_global -= outside_global;
    auto distance = inside_global.two_norm();

    // face volume for integration
    auto facegeo = ig.geometry();
    auto face_volume = facegeo.volume();
    
    // current solution
    auto u_inside = x_i(lfsu_i,0);
    auto u_outside = x_o(lfsu_o,0);

    // bathymmetry
    auto b_inside = Dune::PDELab::Backend::native(bathymmetryvector)[indexset.index(cell_inside)][0];
    auto b_outside = Dune::PDELab::Backend::native(bathymmetryvector)[indexset.index(cell_outside)][0];

    // gradient computed from current solution
    auto vn = (u_inside-u_outside)/distance;

    // permeability
    auto k_inside = model.k(cell_inside,cell_center_local);
    auto k_outside = model.k(cell_outside,cell_center_local);
    auto k = 2.0*(k_inside+1e-15)*(k_outside+1e-15)/(k_inside+k_outside+2e-15);
    
    // nonlinearity
    auto phi = model.phi(b_inside,b_outside,u_inside,u_outside,vn);
    auto dphi_du_inside = model.dphi_du_inside(b_inside,b_outside,u_inside,u_outside,vn);
    auto dphi_du_outside = model.dphi_du_outside(b_inside,b_outside,u_inside,u_outside,vn);
    auto dphi_dvn = model.dphi_dvn(b_inside,b_outside,u_inside,u_outside,vn);
    
    // contribution to jacobian entries
    auto dB_duinside = ((dphi_du_inside+dphi_dvn/distance)*vn + phi/distance)*k*face_volume;
    auto dB_duoutside = ((dphi_du_outside-dphi_dvn/distance)*vn - phi/distance)*k*face_volume;
    mat_ii.accumulate(lfsv_i,0,lfsv_i,0, dB_duinside );
    mat_io.accumulate(lfsv_i,0,lfsv_o,0, dB_duoutside);
    mat_oi.accumulate(lfsv_o,0,lfsv_i,0,-dB_duinside );
    mat_oo.accumulate(lfsv_o,0,lfsv_o,0,-dB_duoutside);
  }

  //! residual contribution of boundary integral (Dirichlet condition)
  // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_boundary (const IG& ig,
                       const LFSU& lfsu_i, const X& x_i, const LFSV& lfsv_i,
		       R& r_i) const
  {
    // check for Dirichlet boundary condition
    auto facegeo = ig.geometry();
    auto facecenterlocal =
      referenceElement(facegeo).position(0,0);
    bool isdirichlet = model.b(ig.intersection(),facecenterlocal);
    if (!isdirichlet) return;

    // inside and outside cells
    auto cell_inside = ig.inside();

    // inside and outside cell geometries
    auto insidegeo = cell_inside.geometry();
    auto cell_center_local = referenceElement(insidegeo).position(0,0);

    // cell centers in global coordinates
    auto inside_global = insidegeo.center();
    auto outside_global = facegeo.center();

    // distance between the two cell centers
    inside_global -= outside_global;
    auto distance = inside_global.two_norm();

    // face volume for integration
    auto face_volume = facegeo.volume();

    // solution
    auto u_inside = x_i(lfsu_i,0);
    auto facecenterinelement=ig.geometryInInside().center();
    auto u_outside = model.g(ig.inside(),facecenterinelement);

    // gradient
    auto vn = (u_inside-u_outside)/distance;

    // permeability
    auto k = model.k(cell_inside,cell_center_local);
    
    // contribution to residual on inside and outside elements
    auto flux = k*vn*face_volume;
    r_i.accumulate(lfsv_i,0, flux);
  }

    //! Jacobian contribution from boundary integral (Dirichlet condition)
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename M>
  void jacobian_boundary (const IG& ig,
                          const LFSU& lfsu_i, const X& x_i,
                          const LFSV& lfsv_i, M& mat_ii) const
  {
    // check for Dirichlet boundary condition
    auto facegeo = ig.geometry();
    auto facecenterlocal =
      referenceElement(facegeo).position(0,0);
    bool isdirichlet = model.b(ig.intersection(),facecenterlocal);
    if (!isdirichlet) return;

    // inside and outside cells
    auto cell_inside = ig.inside();

    // inside and outside cell geometries
    auto insidegeo = cell_inside.geometry();
    auto cell_center_local = referenceElement(insidegeo).position(0,0);

    // cell centers in global coordinates
    auto inside_global = insidegeo.center();
    auto outside_global = facegeo.center();

    // distance between the two cell centers
    inside_global -= outside_global;
    auto distance = inside_global.two_norm();

    // face volume for integration
    auto face_volume = facegeo.volume();

    // solution
    auto u_inside = x_i(lfsu_i,0);
    auto facecenterinelement=ig.geometryInInside().center();
    auto u_outside = model.g(ig.inside(),facecenterinelement);

    // gradient
    auto vn = (u_inside-u_outside)/distance;

    // permeability
    auto k = model.k(cell_inside,cell_center_local);
    
    // contribution to matrix
    mat_ii.accumulate(lfsv_i,0,lfsv_i,0,k*face_volume/distance);
  }

};


/** a local operator for the FV mass operator (L_2 integral)
 *
 * \f{align*}{
 \int_Omega s(x) uv dx
 * \f}
 * \tparam FEM      Type of a finite element map
 */
template<typename Model>
class FVL2
  : public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::FullSkeletonPattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::
  InstationaryLocalOperatorDefaultMethods<double>
{
  Model& model;        // parameter functions
  
public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  enum { doPatternSkeleton = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  //! constructor stores a reference to the parameter object
  FVL2 (Model& model_)
    : model(model_)
  {}

  // set the current time in parameter class
  void setTime (double t)
  {
    model.setTime(t);
  }

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                     const LFSV& lfsv, R& r) const
  {
    // select quadrature rule
    auto geo = eg.geometry();
    auto cellcenterlocal = referenceElement(geo).position(0,0);

    // accumulate residual
    auto s = model.storage(eg.entity(),cellcenterlocal);
    r.accumulate(lfsv,0,s*x(lfsu,0)*geo.volume());
  }

  //! jacobian contribution of volume term
  template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                        M& mat) const
  {
    // types & dimension
    typedef decltype(makeZeroBasisFieldValue(lfsu)) RF;

    // geometry
    auto geo = eg.geometry();
    auto volume = geo.volume();
    auto cellcenterlocal = referenceElement(geo).position(0,0);
    auto s = model.storage(eg.entity(),cellcenterlocal);

    mat.accumulate(lfsv,0,lfsu,0,s*volume);
  }

  //! apply local jacobian of the volume term -> nonlinear variant
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const X& z, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,z,lfsv,r);
  }

  //! apply local jacobian of the volume term -> linear variant
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,x,lfsv,r);
  }

};


// Experimental L2 operator with storage term and piecewise constants
template<typename Model>
class StorageL2 :
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
  const Model& model;
  
public:
  // Pattern assembly flags
  enum { doPatternVolume = true };

  // Residual assembly flags
  enum { doAlphaVolume = true };

  StorageL2 (const Model& model_, int intorderadd=0)
    : model(model_)
    , _intorderadd(intorderadd)
  {}

  // set the current time in parameter class
  void setTime (double time_) const
  {
    model.setTime(time_);
  }

  // prepare step in parameter class
  void preStep (double time_, double dt_, int stages) const
  {
    model.preStep(time_,dt_,stages);
  }

  // Volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // Switches between local and global interface
    using FESwitch = Dune::FiniteElementInterfaceSwitch<
      typename LFSU::Traits::FiniteElementType>;
    using BasisSwitch = Dune::BasisInterfaceSwitch<
      typename FESwitch::Basis>;

    // Define types
    using RF = typename BasisSwitch::RangeField;
    using RangeType = typename BasisSwitch::Range;
    using size_type = typename LFSU::Traits::SizeType;

    // Get geometry
    auto geo = eg.geometry();
    auto cellcenterlocal = referenceElement(geo).position(0,0);
    auto s = model.storage(eg.entity(),cellcenterlocal);

    // Initialize vectors outside for loop
    std::vector<RangeType> phi(lfsu.size());

    // determine integration order
    auto intorder = 2*FESwitch::basis(lfsu.finiteElement()).order() + _intorderadd;

    // Loop over quadrature points
    for (const auto& qp : quadratureRule(geo,intorder))
      {
        // Evaluate basis functions
        FESwitch::basis(lfsu.finiteElement()).evaluateFunction(qp.position(),phi);

        // Evaluate u
        RF u=0.0;
        for (size_type i=0; i<lfsu.size(); i++)
          u += RF(x(lfsu,i)*phi[i]);

        // u*phi_i
        auto factor = s * qp.weight() * geo.integrationElement(qp.position());
        for (size_type i=0; i<lfsu.size(); i++)
          r.accumulate(lfsv,i, u*phi[i]*factor);
      }
  }

  // apply jacobian of volume term
  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
  void jacobian_apply_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
  {
    alpha_volume(eg,lfsu,x,lfsv,y);
  }

  // Jacobian of volume term
  template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & mat) const
  {
    // Switches between local and global interface
    using FESwitch = Dune::FiniteElementInterfaceSwitch<
      typename LFSU::Traits::FiniteElementType>;
    using BasisSwitch = Dune::BasisInterfaceSwitch<
      typename FESwitch::Basis>;

    // Define types
    using RangeType = typename BasisSwitch::Range;
    using size_type = typename LFSU::Traits::SizeType;

    // Get geometry
    auto geo = eg.geometry();
    auto cellcenterlocal = referenceElement(geo).position(0,0);
    auto s = model.storage(eg.entity(),cellcenterlocal);

    // Inititialize vectors outside for loop
    std::vector<RangeType> phi(lfsu.size());

    // determine integration order
    auto intorder = 2*FESwitch::basis(lfsu.finiteElement()).order() + _intorderadd;

    // Loop over quadrature points
    for (const auto& qp : quadratureRule(geo,intorder))
      {
        // Evaluate basis functions
        FESwitch::basis(lfsu.finiteElement()).evaluateFunction(qp.position(),phi);

        // Integrate phi_j*phi_i
        auto factor = s * qp.weight() * geo.integrationElement(qp.position());
        for (size_type j=0; j<lfsu.size(); j++)
          for (size_type i=0; i<lfsu.size(); i++)
            mat.accumulate(lfsv,i,lfsu,j, phi[j]*phi[i]*factor);
      }
  }

private:
  int _intorderadd;
};


#endif
