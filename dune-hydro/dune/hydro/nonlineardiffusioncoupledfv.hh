// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_HYDRO_LOCALOPERATOR_NONLINEARDIFFUSIONCOUPLEDFV_HH
#define DUNE_HYDRO_LOCALOPERATOR_NONLINEARDIFFUSIONCOUPLEDFV_HH

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


template<typename T>
T smoothmax0 (T x, T tau)
{
  if (x<=0.0) return 0.0;
  if (x>=tau) return x;
  return tau*(x/tau)*(x/tau)*(2.0-x/tau);
}

template<typename T>
T smoothmin0 (T x, T tau)
{
  if (x<=-tau) return x;
  if (x>=0.0) return 0.0;
  return -tau*(x/tau)*(x/tau)*(2.0+x/tau);
}

/** a local operator for spatial part of two coupled nonlinear diffusion equations with cell-centered FV
 *  The initial version of this operator has been implemented by Johanna Biereder during her
 *  Bachelor thesis in 2021.
 *
 */
template<typename Model, typename GFS>
class NonlinearDiffusionCoupledFV :
  public Dune::PDELab::NumericalJacobianVolume<NonlinearDiffusionCoupledFV<Model,GFS> >,
  public Dune::PDELab::NumericalJacobianSkeleton<NonlinearDiffusionCoupledFV<Model,GFS> >,
  public Dune::PDELab::NumericalJacobianBoundary<NonlinearDiffusionCoupledFV<Model,GFS> >,
  public Dune::PDELab::NumericalJacobianApplyVolume<NonlinearDiffusionCoupledFV<Model,GFS> >,
  public Dune::PDELab::NumericalJacobianApplySkeleton<NonlinearDiffusionCoupledFV<Model,GFS> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<NonlinearDiffusionCoupledFV<Model,GFS> >,
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
  Z bathymmetryvector_surface;
  Z bathymmetryvector_groundwater;

public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  enum { doPatternSkeleton = true };

  // residual assembly flags
  enum { doLambdaVolume = true };
  enum { doLambdaBoundary = true };
  enum { doAlphaSkeleton  = true };
  enum { doAlphaBoundary  = true };
  enum { doAlphaVolume  = true };

  //! constructor stores a reference to the parameter object
  NonlinearDiffusionCoupledFV (Model& model_, const GFS& gfs_)
    : Dune::PDELab::NumericalJacobianVolume<NonlinearDiffusionCoupledFV<Model,GFS>>(1e-10),
    Dune::PDELab::NumericalJacobianSkeleton<NonlinearDiffusionCoupledFV<Model,GFS> >(1e-10),
    Dune::PDELab::NumericalJacobianBoundary<NonlinearDiffusionCoupledFV<Model,GFS> >(1e-10),
    model(model_), pgfs(stackobject_to_shared_ptr(gfs_)), gv(pgfs->gridView()), 
      indexset(gv.indexSet()), bathymmetryvector_surface(*pgfs), bathymmetryvector_groundwater(*pgfs)
  {
    auto bathymmetrylambda_surface = [&](const auto& e, const auto& x)
      {return model.bathymmetry_surface(e,x);};
    auto bathymmetrylambda_groundwater = [&](const auto& e, const auto& x)
      {return model.bathymmetry_groundwater(e,x);};
    auto bathymmetrygf_surface = Dune::PDELab::
      makeInstationaryGridFunctionFromCallable(gv,bathymmetrylambda_surface,model);
    auto bathymmetrygf_groundwater = Dune::PDELab::
      makeInstationaryGridFunctionFromCallable(gv,bathymmetrylambda_groundwater,model);
    Dune::PDELab::interpolate(bathymmetrygf_surface,*pgfs,bathymmetryvector_surface);
    Dune::PDELab::interpolate(bathymmetrygf_groundwater,*pgfs,bathymmetryvector_groundwater);
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
    using namespace Dune::Indices;  
    auto lfsv_surface = lfsv.child(_0);
    auto lfsv_groundwater = lfsv.child(_1);

    // center of reference element
    auto cellgeo = eg.geometry();
    auto cellcenterlocal =
      referenceElement(cellgeo).position(0,0);

    // accumulate residual
    auto f_surface = model.f_surface(eg.entity(),cellcenterlocal);
    auto f_groundwater = model.f_groundwater(eg.entity(),cellcenterlocal);
    r.accumulate(lfsv_surface,0,-f_surface*cellgeo.volume());
    r.accumulate(lfsv_groundwater,0,-f_groundwater*cellgeo.volume());
  }

  //! residual contribution of rhs boundary integral (Neumann boundary condition)
  template<typename IG, typename LFSV, typename R>
  void lambda_boundary (const IG& ig, const LFSV& lfsv_i,   
                        R& r_i) const
  {
    using namespace Dune::Indices;  
    auto lfsv_i_surface = lfsv_i.child(_0);
    auto lfsv_i_groundwater = lfsv_i.child(_1);

    // face volume for integration
    auto facegeo = ig.geometry();
    auto facecenterlocal = referenceElement(facegeo).position(0,0);

    // surface Neumann condition
    bool isdirichlet_surface = model.b_surface(ig.intersection(),facecenterlocal);
    if (!isdirichlet_surface) {
      // contribution to residual from Neumann boundary
      auto j_groundwater = model.j_groundwater(ig.intersection(),facecenterlocal);
      r_i.accumulate(lfsv_i_groundwater,0,j_groundwater*facegeo.volume());
    }

    // groundwater Neumann condition
    bool isdirichlet_groundwater = model.b_groundwater(ig.intersection(),facecenterlocal);
    if (!isdirichlet_groundwater) {
      // contribution to residual from Neumann boundary
      auto j_surface = model.j_surface(ig.intersection(),facecenterlocal);
      r_i.accumulate(lfsv_i_surface,0,j_surface*facegeo.volume());
    }
  }

  //! residual contribution of volume integral (reaction term coupling the two equations)
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                     const LFSV& lfsv, R& r) const
  {
    using namespace Dune::Indices;  
    auto lfsu_surface = lfsu.child(_0);
    auto lfsu_groundwater = lfsu.child(_1);

    auto lfsv_surface = lfsv.child(_0);
    auto lfsv_groundwater = lfsv.child(_1);

    // center of reference element
    auto cellgeo = eg.geometry();
    auto cellcenterlocal =
      referenceElement(cellgeo).position(0,0);

    // get cell values
    auto u_surface = x(lfsu_surface,0);
    auto u_groundwater = x(lfsu_groundwater,0);
    auto bathymmetry_surface = model.bathymmetry_surface(eg.entity(),cellcenterlocal);

    // evaluate reaction term
    auto q = model.q(eg.entity(),cellcenterlocal,u_surface,u_groundwater,bathymmetry_surface);

    // and accumulate
    r.accumulate(lfsv_surface,0,q*eg.geometry().volume());        
    r.accumulate(lfsv_groundwater,0,-q*eg.geometry().volume());    
  }

  
  //! residual contribution from skeleton terms
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_skeleton (const IG& ig,
         const LFSU& lfsu_i, const X& x_i, const LFSV& lfsv_i,
         const LFSU& lfsu_o, const X& x_o, const LFSV& lfsv_o,
         R& r_i, R& r_o) const
  {
    using namespace Dune::Indices;
    auto lfsu_i_surface = lfsu_i.child(_0);
    auto lfsu_i_groundwater = lfsu_i.child(_1);
    auto lfsu_o_surface = lfsu_o.child(_0);
    auto lfsu_o_groundwater = lfsu_o.child(_1);
 
    auto lfsv_i_surface = lfsv_i.child(_0);
    auto lfsv_i_groundwater = lfsv_i.child(_1);
    auto lfsv_o_surface = lfsv_o.child(_0);
    auto lfsv_o_groundwater = lfsv_o.child(_1);

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
    auto u_inside_surface = x_i(lfsu_i_surface,0);
    auto u_outside_surface = x_o(lfsu_o_surface,0);
    auto u_inside_groundwater = x_i(lfsu_i_groundwater,0);
    auto u_outside_groundwater = x_o(lfsu_o_groundwater,0);

    // bathymmetry
    // auto b_inside_surface = Dune::PDELab::Backend::native(bathymmetryvector_surface)[indexset.index(cell_inside)][0];
    // auto b_outside_surface = Dune::PDELab::Backend::native(bathymmetryvector_surface)[indexset.index(cell_outside)][1];
    // auto b_inside_groundwater = Dune::PDELab::Backend::native(bathymmetryvector_groundwater)[indexset.index(cell_inside)][0];
    // auto b_outside_groundwater = Dune::PDELab::Backend::native(bathymmetryvector_groundwater)[indexset.index(cell_outside)][1];
    auto b_inside_surface = model.bathymmetry_surface(cell_inside,cell_center_local);
    auto b_outside_surface = model.bathymmetry_surface(cell_outside,cell_center_local);
    auto b_inside_groundwater = model.bathymmetry_groundwater(cell_inside,cell_center_local);
    auto b_outside_groundwater = model.bathymmetry_groundwater(cell_outside,cell_center_local);

    // gradient computed from current solution
    auto vn_surface = (u_inside_surface-u_outside_surface)/distance;
    auto vn_groundwater = (u_inside_groundwater-u_outside_groundwater)/distance;

    // permeability
    auto k_inside_surface = model.k_surface(cell_inside,cell_center_local);
    auto k_outside_surface = model.k_surface(cell_outside,cell_center_local);
    auto k_surface = 2.0*(k_inside_surface+1e-15)*(k_outside_surface+1e-15)/(k_inside_surface+k_outside_surface+2e-15);

    auto k_inside_groundwater = model.k_groundwater(cell_inside,cell_center_local);
    auto k_outside_groundwater = model.k_groundwater(cell_outside,cell_center_local);
    auto k_groundwater = 2.0*(k_inside_groundwater+1e-15)*(k_outside_groundwater+1e-15)/(k_inside_groundwater+k_outside_groundwater+2e-15);
  
    // nonlinearity
    auto phi_surface = model.phi_surface(b_inside_surface,b_outside_surface,u_inside_surface,u_outside_surface,vn_surface);
    auto phi_groundwater = model.phi_groundwater(b_inside_groundwater,b_outside_groundwater,u_inside_groundwater,u_outside_groundwater,vn_groundwater);
    
    // contribution to residual on inside and outside elements
    auto flux_surface = k_surface*phi_surface*vn_surface*face_volume;
    auto flux_groundwater = k_groundwater*phi_groundwater*vn_groundwater*face_volume;

    // std::cout << "inside cell #" << indexset.index(cell_inside) << " outside cell #" << indexset.index(cell_outside) << std::endl;
    // std::cout << "  surface: b_i=" << b_inside_surface << " b_o=" <<  b_outside_surface << " u_i=" << u_inside_surface << " u_o=" << u_outside_surface
    //           << " vn=" << vn_surface << " phi=" << phi_surface << std::endl;
    // std::cout << "  groundw: b_i=" << b_inside_groundwater << " b_o=" <<  b_outside_groundwater << " u_i=" << u_inside_groundwater << " u_o=" << u_outside_groundwater
    //           << " vn=" << vn_groundwater << " phi=" << phi_groundwater << std::endl;
    
    r_i.accumulate(lfsv_i_surface,0, flux_surface);
    r_o.accumulate(lfsv_o_surface,0,-flux_surface);
    r_i.accumulate(lfsv_i_groundwater,0, flux_groundwater);
    r_o.accumulate(lfsv_o_groundwater,0,-flux_groundwater);
  }

  //! Jacobian contribution from skeleton terms - not working, error in accumulate
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename M>
  void xxx_jacobian_skeleton (const IG& ig,
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

    // bathymmetry    
    auto b_inside_surface = model.bathymmetry_surface(cell_inside,cell_center_local);
    auto b_outside_surface = model.bathymmetry_surface(cell_outside,cell_center_local);

    auto b_inside_groundwater = model.bathymmetry_groundwater(cell_inside,cell_center_local);
    auto b_outside_groundwater = model.bathymmetry_groundwater(cell_outside,cell_center_local);

    // access child spaces ansatz/test inside/outside surface/groundwater = 8 local function spaces
    using namespace Dune::Indices;   
    auto lfsu_i_surface = lfsu_i.child(_0);
    auto lfsu_i_groundwater = lfsu_i.child(_1);
    auto lfsu_o_surface = lfsu_o.child(_0);
    auto lfsu_o_groundwater = lfsu_o.child(_1);
    auto lfsv_i_surface = lfsv_i.child(_0);
    auto lfsv_i_groundwater = lfsv_i.child(_1);
    auto lfsv_o_surface = lfsv_o.child(_0);
    auto lfsv_o_groundwater = lfsv_o.child(_1);

    // surface
    auto u_inside_surface = x_i(lfsu_i_surface,0);
    auto u_outside_surface = x_o(lfsu_o_surface,0);

    // gradient computed from current solution
    auto vn_surface = (u_inside_surface-u_outside_surface)/distance;

    // permeability
    auto k_inside_surface = model.k_surface(cell_inside,cell_center_local);
    auto k_outside_surface = model.k_surface(cell_outside,cell_center_local);
    auto k_surface = 2.0*(k_inside_surface+1e-15)*(k_outside_surface+1e-15)/(k_inside_surface+k_outside_surface+2e-15);
    
    // nonlinearity
    auto phi_surface = model.phi_surface(b_inside_surface,b_outside_surface,u_inside_surface,u_outside_surface,vn_surface);
    auto dphi_du_inside_surface = model.dphi_du_inside_surface(b_inside_surface,b_outside_surface,u_inside_surface,u_outside_surface,vn_surface);
    auto dphi_du_outside_surface = model.dphi_du_outside_surface(b_inside_surface,b_outside_surface,u_inside_surface,u_outside_surface,vn_surface);
    auto dphi_dvn_surface = model.dphi_dvn_surface(b_inside_surface,b_outside_surface,u_inside_surface,u_outside_surface,vn_surface);
    
    // contribution to jacobian entries
    auto dB_duinside_surface = ((dphi_du_inside_surface+dphi_dvn_surface/distance)*vn_surface + phi_surface/distance)*k_surface*face_volume;
    auto dB_duoutside_surface = ((dphi_du_outside_surface-dphi_dvn_surface/distance)*vn_surface - phi_surface/distance)*k_surface*face_volume;
    mat_ii.accumulate(lfsv_i_surface,0,lfsv_i_surface,0, dB_duinside_surface );
    mat_io.accumulate(lfsv_i_surface,0,lfsv_o_surface,0, dB_duoutside_surface);
    mat_oi.accumulate(lfsv_o_surface,0,lfsv_i_surface,0,-dB_duinside_surface );
    mat_oo.accumulate(lfsv_o_surface,0,lfsv_o_surface,0,-dB_duoutside_surface);

    // groundwater
    auto u_inside_groundwater = x_i(lfsu_i_groundwater,0);
    auto u_outside_groundwater = x_o(lfsu_o_groundwater,0);

    // gradient computed from current solution
    auto vn_groundwater = (u_inside_groundwater-u_outside_groundwater)/distance;

    // permeability
    auto k_inside_groundwater = model.k_groundwater(cell_inside,cell_center_local);
    auto k_outside_groundwater = model.k_groundwater(cell_outside,cell_center_local);
    auto k_groundwater = 2.0*(k_inside_groundwater+1e-15)*(k_outside_groundwater+1e-15)/(k_inside_groundwater+k_outside_groundwater+2e-15);
    
    // nonlinearity
    auto phi_groundwater = model.phi_groundwater(b_inside_groundwater,b_outside_groundwater,u_inside_groundwater,u_outside_groundwater,vn_groundwater);
    auto dphi_du_inside_groundwater = model.dphi_du_inside_groundwater(b_inside_groundwater,b_outside_groundwater,u_inside_groundwater,u_outside_groundwater,vn_groundwater);
    auto dphi_du_outside_groundwater = model.dphi_du_outside_groundwater(b_inside_groundwater,b_outside_groundwater,u_inside_groundwater,u_outside_groundwater,vn_groundwater);
    auto dphi_dvn_groundwater = model.dphi_dvn_groundwater(b_inside_groundwater,b_outside_groundwater,u_inside_groundwater,u_outside_groundwater,vn_groundwater);
    
    // contribution to jacobian entries
    auto dB_duinside_groundwater = ((dphi_du_inside_groundwater+dphi_dvn_groundwater/distance)*vn_groundwater + phi_groundwater/distance)*k_groundwater*face_volume;
    auto dB_duoutside_groundwater = ((dphi_du_outside_groundwater-dphi_dvn_groundwater/distance)*vn_groundwater - phi_groundwater/distance)*k_groundwater*face_volume;
    mat_ii.accumulate(lfsv_i_groundwater,0,lfsv_i_groundwater,0, dB_duinside_groundwater );
    mat_io.accumulate(lfsv_i_groundwater,0,lfsv_o_groundwater,0, dB_duoutside_groundwater);
    mat_oi.accumulate(lfsv_o_groundwater,0,lfsv_i_groundwater,0,-dB_duinside_groundwater );
    mat_oo.accumulate(lfsv_o_groundwater,0,lfsv_o_groundwater,0,-dB_duoutside_groundwater);
  } 

  //! residual contribution of boundary integral (Dirichlet condition)
  // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
  template<typename IG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_boundary (const IG& ig,
                       const LFSU& lfsu_i, const X& x_i, const LFSV& lfsv_i,
                       R& r_i) const
  {
    using namespace Dune::Indices;   
    auto lfsu_i_surface = lfsu_i.child(_0);
    auto lfsu_i_groundwater = lfsu_i.child(_1);
    auto lfsv_i_surface = lfsv_i.child(_0);
    auto lfsv_i_groundwater = lfsv_i.child(_1);

    // check for Dirichlet boundary condition
    auto facegeo = ig.geometry();
    auto facecenterlocal = referenceElement(facegeo).position(0,0);
    bool isdirichlet_surface = model.b_surface(ig.intersection(),facecenterlocal);
    bool isdirichlet_groundwater = model.b_groundwater(ig.intersection(),facecenterlocal);
    if (isdirichlet_surface==false && isdirichlet_groundwater==false) return;

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

    // For solution
    auto facecenterinelement=ig.geometryInInside().center();

    if (isdirichlet_surface)
      {
        // solution
        auto u_inside_surface = x_i(lfsu_i_surface,0);
        auto u_outside_surface = model.g_surface(ig.inside(),facecenterinelement);

        // gradient
        auto vn_surface = (u_inside_surface-u_outside_surface)/distance;

        // permeability
        auto k_surface = model.k_surface(cell_inside,cell_center_local);
    
        // contribution to residual on inside and outside elements
        auto flux_surface = k_surface*vn_surface*face_volume;
        r_i.accumulate(lfsv_i_surface,0, flux_surface);
      }

    if (isdirichlet_groundwater)
      {
        // solution
        auto u_inside_groundwater = x_i(lfsu_i_groundwater,0);
        auto u_outside_groundwater = model.g_groundwater(ig.inside(),facecenterinelement);

        // gradient
        auto vn_groundwater = (u_inside_groundwater-u_outside_groundwater)/distance;

        // permeability
        auto k_groundwater = model.k_groundwater(cell_inside,cell_center_local);
    
        // contribution to residual on inside and outside elements
        auto flux_groundwater = k_groundwater*vn_groundwater*face_volume;
        r_i.accumulate(lfsv_i_groundwater,0, flux_groundwater);
      }
  }
};


/** a local operator for the FV mass operator (L_2 integral)
 *
 * \f{align*}{
 \int_Omega s(x) uv dx
 * \f}
 * \tparam FEM      Type of a finite element map
 */
template<typename Model, typename GFS>
class FVTemporalCoupled
  : public Dune::PDELab::NumericalJacobianVolume<FVTemporalCoupled<Model,GFS> >,
    public Dune::PDELab::NumericalJacobianApplyVolume<FVTemporalCoupled<Model,GFS> >,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::FullSkeletonPattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::
  InstationaryLocalOperatorDefaultMethods<double>
{
  Model& model;        // parameter functions

  typedef typename GFS::Traits::GridView GV;
  typedef typename GV::IndexSet IndexSet;
  using Z = Dune::PDELab::Backend::Vector<GFS,double>;
  std::shared_ptr<const GFS> pgfs;
  const GV& gv;
  const IndexSet& indexset;
  Z bathymmetryvector_surface;
  Z bathymmetryvector_groundwater;

public:
  // pattern assembly flags
  enum { doPatternVolume = true };
  enum { doPatternSkeleton = true };

  // residual assembly flags
  enum { doAlphaVolume = true };

  //! constructor stores a reference to the parameter object
  FVTemporalCoupled (Model& model_, const GFS& gfs_)
    : Dune::PDELab::NumericalJacobianVolume<FVTemporalCoupled<Model,GFS> >(1e-10),
    model(model_), pgfs(stackobject_to_shared_ptr(gfs_)), gv(pgfs->gridView()), 
      indexset(gv.indexSet()), bathymmetryvector_surface(*pgfs), bathymmetryvector_groundwater(*pgfs)
  {
    auto bathymmetrylambda_surface = [&](const auto& e, const auto& x)
      {return model.bathymmetry_surface(e,x);};
    auto bathymmetrylambda_groundwater = [&](const auto& e, const auto& x)
      {return model.bathymmetry_groundwater(e,x);};
    auto bathymmetrygf_surface = Dune::PDELab::
      makeInstationaryGridFunctionFromCallable(gv,bathymmetrylambda_surface,model);
    auto bathymmetrygf_groundwater = Dune::PDELab::
      makeInstationaryGridFunctionFromCallable(gv,bathymmetrylambda_groundwater,model);
    Dune::PDELab::interpolate(bathymmetrygf_surface,*pgfs,bathymmetryvector_surface);
    Dune::PDELab::interpolate(bathymmetrygf_groundwater,*pgfs,bathymmetryvector_groundwater);
  }

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
    
    // select the two components
    using namespace Dune::Indices;
    auto lfsu_surface = lfsu.child(_0);
    auto lfsu_groundwater = lfsu.child(_1);
    auto lfsv_surface = lfsv.child(_0);
    auto lfsv_groundwater = lfsv.child(_1);

    // select quadrature rule
    auto geo = eg.geometry();
    auto cellcenterlocal = referenceElement(geo).position(0,0);

    // storage coefficient
    auto storage_surface = model.storage_surface(eg.entity(),cellcenterlocal);
    auto storage_groundwater = model.storage_groundwater(eg.entity(),cellcenterlocal);

    // current solution
    auto u_surface = x(lfsu_surface,0);
    auto u_groundwater = x(lfsu_groundwater,0);

    // bathymmetry
    // auto bathymmetry_surface = Dune::PDELab::Backend::native(bathymmetryvector_surface)[indexset.index(eg.entity())][0];
    // auto bathymmetry_groundwater = Dune::PDELab::Backend::native(bathymmetryvector_groundwater)[indexset.index(eg.entity())][0];
    auto bathymmetry_surface = model.bathymmetry_surface(eg.entity(),cellcenterlocal);
    auto bathymmetry_groundwater = model.bathymmetry_groundwater(eg.entity(),cellcenterlocal);

    double value = smoothmin0(u_groundwater-bathymmetry_surface,1.0) + (bathymmetry_surface - bathymmetry_groundwater);

    // accumulate residual
    r.accumulate(lfsv_surface,0,storage_surface*(u_surface-bathymmetry_surface)*geo.volume());
    r.accumulate(lfsv_groundwater,0,storage_groundwater*value*geo.volume());
  }
};

#endif
