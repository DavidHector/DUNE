// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// always include the config file
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
// C++ includes
#include<math.h>
#include<iostream>
#include<vector>
#include<set>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
#include <dune/common/exceptions.hh> // We use exceptions
// dune-istl includes
#include<dune/istl/preconditioners.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/onedgrid.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/alugrid/dgf.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
// dune-pdelab includes
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/boilerplate/pdelab.hh>
// ftworth includes
#include"../dune/ftworth/partitioner.hh"
#include"../dune/ftworth/assemblewrapper.hh"
#include"../dune/ftworth/subdomainutilities.hh"
#include"../dune/ftworth/cdediscretizer.hh"
#include"../dune/ftworth/schwarz_preconditioners.hh"
#include"../dune/ftworth/twolevel_geneo_preconditioner.hh"
// includes from this directory
#include"poisson.hh"
#include"heterogeneous_problem.hh"


int main(int argc, char** argv)
{
  try{
    // Maybe initialize MPI
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    std::cout << "Hello World! This is dune-ftworth." << std::endl;
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()
               <<" processes!"<<std::endl;

    // Read parameters from ini file
    Dune::ParameterTree ptree;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("DG_diagnostics.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);

#if (HAVE_UG) // check if we have the UG grid manager

    //===================
    // grid generation and partitioning
    //===================
    
    // load a grid file and refine
    std::string mshfile = ptree.get("grid.filename","circle.msh");
    using Grid = Dune::UGGrid<2>;  
    std::shared_ptr<Grid> pgrid;
    if (helper.rank()==0)
      // the third argument must be false, otherwise it does not work ... not documented!
      pgrid = std::shared_ptr<Grid>(Dune::GmshReader<Grid>::read(mshfile,true,false));
    
    // refine grid
    int refine = ptree.get("grid.refine",(int)0);
    pgrid->globalRefine(ptree.get("grid.refine",refine));

    // // create structured grid
    // typedef Dune::YaspGrid<2> GM;
    // typedef Dune::PDELab::StructuredGrid<GM> Grid;
    // Grid pgrid(Dune::GeometryType::cube,16,1);
    // pgrid->globalRefine(2);

    // extract the leaf grid view
    auto gv = pgrid->leafGridView();
    using GV = decltype(gv); // and get the type
    const int dim = GV::dimension;

    // partition the elements using ParMetis
    unsigned subdomains = ptree.get("geneo.subdomains",(unsigned)4);
    auto partition = parmetis_partitioning(gv,helper,subdomains);

    // generate overlap: need to decide on size of overlap as well as extension method
    unsigned overlap = ptree.get("geneo.overlap",(unsigned)1);
    std::string extensionmethod = ptree.get("geneo.extensionmethod","vertex");
    auto overlappingsubdomains = grow_subdomains(gv,subdomains,partition,overlap,extensionmethod);

    // test overlap generation
    auto& indexset = gv.indexSet(); // to attach data to elements
    std::cout << indexset.size(0) << " elements" << std::endl;
    std::cout << indexset.size(dim) << " vertices" << std::endl;
    std::cout << subdomains << " subdomains" << std::endl;
    int vp = ptree.get("geneo.view",(int)0);
    std::vector<unsigned> adomain(indexset.size(0));
    std::vector<unsigned> k0(indexset.size(0));
    size_t k_0 = 0;
    unsigned count_elements=0;
    for (const auto& e : elements(gv))
      {
        adomain[indexset.index(e)] = 0;
        if (overlappingsubdomains[indexset.index(e)].count(vp)>0) adomain[indexset.index(e)] = 1;
        if (partition[indexset.index(e)]==vp) adomain[indexset.index(e)] = 2;
        k0[indexset.index(e)] = overlappingsubdomains[indexset.index(e)].size();
        k_0 = std::max(k_0,overlappingsubdomains[indexset.index(e)].size());
        count_elements += overlappingsubdomains[indexset.index(e)].size();
      }
    std::cout << k_0 << " k_0" << std::endl;
    std::cout << indexset.size(0) << " elements, " << count_elements << " in all subdomains, ratio=" << ((double)count_elements)/indexset.size(0) << std::endl;

    //===================
    // set up a discretization scheme
    //===================

    // solve a Poisson problem
    using DF = Grid::ctype;
    using RF = double;

    // set up the solver class
    using Problem = PoissonProblem<GV,RF>;
    // using Problem = HeterogeneousProblem<GV,RF>;    
    Problem problem;
    using Discretization = DGConvectionDiffusionProblemSimplex<GV,Problem,1>;
    Discretization disc(gv,problem,ptree);

    //===================
    // assemble the global problem
    //===================

    // set up linear system in PDELab components
    using VEC = typename Discretization::VEC;
    using MAT = typename Discretization::MAT;

    VEC x(*disc.gfs);
    Dune::PDELab::interpolate(*disc.g,*disc.gfs,x);
    MAT A(*disc.go);
    A = 0.0;
    disc.go->jacobian(x,A);
    VEC r(*disc.gfs);
    r = 0.0;
    disc.go->residual(x,r);

    // extract ISTL objects
    using ISTLV = Dune::PDELab::Backend::Native<VEC>;
    ISTLV& istlx = Dune::PDELab::Backend::native(x);
    ISTLV& istlr = Dune::PDELab::Backend::native(r);
    using ISTLM = Dune::PDELab::Backend::Native<MAT>;
    ISTLM& istlA = Dune::PDELab::Backend::native(A);
    std::shared_ptr<ISTLM> pistlA = stackobject_to_shared_ptr(Dune::PDELab::Backend::native(istlA));

    //===================
    // build subdomain information
    //==================

    std::cout << "building subdomain information" << std::endl;
    using GFS = typename Discretization::GFS;
    using FLSDI = FineLevelSubdomainInformation<GFS>;
    FLSDI flsdi(*disc.gfs,subdomains,overlappingsubdomains);
    unsigned summed_dofs=0;
    for (unsigned i=0; i<flsdi.global_to_subdomains->size(); i++)
      summed_dofs += (*flsdi.global_to_subdomains)[i].size();
    std::cout << flsdi.global_to_subdomains->size() << " dof blocks, " << summed_dofs << " dof blocks in all subdomains, ratio=" << ((double)summed_dofs)/flsdi.global_to_subdomains->size() << std::endl;

    unsigned dofspersubdomain=0;
    for (unsigned p=0; p<subdomains; ++p)
      if ((*flsdi.local_to_global)[p].size()>dofspersubdomain) dofspersubdomain = (*flsdi.local_to_global)[p].size();
    std::cout << "maximum number of dofs (blocks) per subdomain " << dofspersubdomain << std::endl;

    //===================
    // set up subdomain matrices
    //==================

    std::cout << "setting up subdomain matrices" << std::endl;
    auto matrices = std::shared_ptr<std::vector<std::shared_ptr<ISTLM>>>(new std::vector<std::shared_ptr<ISTLM>>(set_up_subdomain_matrices(istlA,*flsdi.local_to_global,*flsdi.global_to_subdomains,*flsdi.global_to_local)));

    // check subdomain matrices
    // for (unsigned p=0; p<subdomains; ++p)
    //   {
    //     ISTLM& myA = *matrices[p];
    //     std::cout << "matrix in subdomain " << p << " size " << myA.N() << " times " << myA.M() << std::endl;
    //     for (size_t i=0; i<myA.N(); i++)
    //     {
    //       auto iglobal = flsdi.local_to_global[p][i];
    //       std::cout << "row [" << i << "|" << iglobal << "] =";
    //       auto cIt = myA[i].begin();
    //       auto cEndIt = myA[i].end();
    //       for (; cIt!=cEndIt; ++cIt)
    //         {
    //           auto j = cIt.index();
    //           auto jglobal = flsdi.local_to_global[p][j];
    //           std::cout << " (" << j << "|" << jglobal << ")";
    //         }
    //       std::cout << std::endl;
    //     }
    //   }

    //===================
    // assemble subdomain matrices
    //==================

    std::cout << "assembling subdomain matrices" << std::endl;
    // now we may wrap the local operator
    using WLOP = SubdomainAssemblerWrapper<typename Discretization::LOP,ISTLM,FLSDI>;
    WLOP wlop(*disc.lop,overlappingsubdomains,flsdi,matrices);

    // and set up a new grid operator
    using WGOP = Dune::PDELab::GridOperator<typename Discretization::GFS,typename Discretization::GFS,WLOP,
                                            typename Discretization::MBE,typename Discretization::RF,typename Discretization::RF,typename Discretization::RF,
                                            typename Discretization::CC,typename Discretization::CC>;
    WGOP wgop(*disc.gfs,*disc.cc,*disc.gfs,*disc.cc,wlop,*disc.mbe);
    A = 0.0;
    for (unsigned p=0; p<subdomains; ++p)
      *(*matrices)[p] = 0.0; // clear subdomain matrices
    wgop.jacobian(x,A);

    // still need to set constraints in subdomain matrices
    // we do the same trick: make a coefficient vector indicating where the constraints are
    std::cout << "assembling dirichlet constraints" << std::endl;
    VEC dirichlet(*disc.gfs);
    dirichlet = 1.0;
    set_constrained_dofs(*disc.cc,0.0,dirichlet); // dirichlet dofs are zero, others are 1, so we can use it as a mask!
    std::shared_ptr<ISTLV> pistldirichlet = stackobject_to_shared_ptr(Dune::PDELab::Backend::native(dirichlet));
    
    auto floating = std::shared_ptr<std::vector<bool>>(new std::vector<bool>(assemble_dirichlet_constraints(*flsdi.local_to_global,Dune::PDELab::Backend::native(dirichlet),*matrices)));
    int floatingsubdomains=0;
    for (unsigned p=0; p<subdomains; ++p)
      if ((*floating)[p])
        floatingsubdomains++;
    std::cout << floatingsubdomains << " floating subdomains" << std::endl;        
    
    // now we know where the Dirichlet constraints for subdomain boundaries are ...
    auto boundary_dofs = std::shared_ptr<std::vector<std::vector<bool>>>(new std::vector<std::vector<bool>>(get_subdomain_boundary_dofs(istlA,*flsdi.local_to_global,*flsdi.global_to_subdomains)));
    std::cout << "calculated boundary dofs" << std::endl;

    // // compare subdomain matrices and global matrix
    // std::cout << "comparing matrices" << std::endl;
    // for (unsigned p=0; p<subdomains; ++p)
    //   {
    //     // check dimensions
    //     if ( (*matrices)[p]->N()!=istlA.N() || (*matrices)[p]->M()!=istlA.M() )
    //       {
    //         std::cout << p << ": matrix size differs " << (*matrices)[p]->N() << "x" << (*matrices)[p]->M() << " vs " << istlA.N() << "x" << istlA.M() << std::endl;
    //         continue;
    //       }
    //     // compare matrix
    //     ISTLM& Asubdomain = *(*matrices)[p]; // modify the matrix in place
    //     auto blocksize = flsdi.blocksize;
    //     unsigned differ = 0;
    //     unsigned count = 0;
    //     for (size_t i=0; i<istlA.N(); i++)
    //       {
    //         auto isubdomain = (*flsdi.global_to_local)[i][p];
    //         auto cIt = istlA[i].begin();
    //         auto cEndIt = istlA[i].end();
    //         for (; cIt!=cEndIt; ++cIt)
    //           {
    //             auto j = cIt.index();
    //             auto jsubdomain = (*flsdi.global_to_local)[j][p];
    //             for (int compi=0; compi<blocksize; compi++)
    //               for (int compj=0; compj<blocksize; compj++)
    //                 {
    //                   typename ISTLM::block_type A1 = *cIt;
    //                   typename ISTLM::block_type A2 = Asubdomain[isubdomain][jsubdomain];
    //                   count++;
    //                   if ( std::abs(A1[compi][compj]-A2[compi][compj]) > 1e-12 )
    //                     {
    //                       std::cout << "!" << p << ": " << i << "," << j << "," << compi << "," << compj << " = " << A1[compi][compj] << " : " << A2[compi][compj] << std::endl;
    //                       differ++;
    //                     }
    //                   // else
    //                   //   {
    //                   //     std::cout << "=" << p << ": " << i << "," << j << "," << compi << "," << compj << " = " << A1[compi][compj] << " : " << A2[compi][compj] << std::endl;
    //                   //   }
    //                 }
    //           }
    //       }
    //     if (differ)
    //       std::cout << p << ": " << differ << " different matrix entries out of " << count << std::endl;
    //     else
    //       std::cout << p << ": " << " all " << count << " matrix entries coincide" << std::endl;
    //   }

    
    //===================
    // now solve the problem
    //==================

    // display partition of unity
    std::string pum = ptree.get("geneo.pum","standard");
    std::vector<std::vector<RF>> partition_of_unity;
    if (pum=="distance")
      partition_of_unity = get_partition_of_unity_distance(matrices,flsdi.local_to_global,flsdi.global_to_local,boundary_dofs);
    else
      partition_of_unity = get_partition_of_unity_standard<RF>(flsdi.local_to_global,flsdi.global_to_local,boundary_dofs);
    VEC puvec(*disc.gfs);
    puvec = 0.0;
    ISTLV& istlpuvec = Dune::PDELab::Backend::native(puvec);
    for (unsigned i=0; i<(*flsdi.local_to_global)[vp].size(); i++)
      istlpuvec[(*flsdi.local_to_global)[vp][i]] = partition_of_unity[vp][i];
    std::cout << "made PU" << std::endl;

    // make coarse grid vectors for two level method
    const int nvec = ptree.get<int>("geneo.number_of_basisvectors",3);
    std::cout << nvec << "vectors" << std::endl;
    std::vector<std::shared_ptr<VEC>> coarse_vectors(nvec);
    for (int i=0; i<nvec; i++) coarse_vectors[i] = std::shared_ptr<VEC>(new VEC(*disc.gfs,0.0));
    std::vector<std::shared_ptr<ISTLV>> istl_coarse_vectors(nvec);
    for (int i=0; i<nvec; i++) istl_coarse_vectors[i] = stackobject_to_shared_ptr(Dune::PDELab::Backend::native(*coarse_vectors[i]));
      
    auto constlambda = [](const auto& x){ return 1.0;};
    auto constfunc = Dune::PDELab::makeGridFunctionFromCallable(gv,constlambda);
    Dune::PDELab::interpolate(constfunc,*disc.gfs,*coarse_vectors[0]);

    if (nvec>=2)
    {
      auto xlambda = [](const auto& x){ return x[0];};
      auto xfunc = Dune::PDELab::makeGridFunctionFromCallable(gv,xlambda);
      Dune::PDELab::interpolate(xfunc,*disc.gfs,*coarse_vectors[1]);
    }
    if (nvec>=3)
    {
      auto ylambda = [](const auto& x){ return x[1];};
      auto yfunc = Dune::PDELab::makeGridFunctionFromCallable(gv,ylambda);
      Dune::PDELab::interpolate(yfunc,*disc.gfs,*coarse_vectors[2]);
    }
    if (nvec>=4 && dim==3)
      {
        auto zlambda = [](const auto& x){ return x[2];};
        auto zfunc = Dune::PDELab::makeGridFunctionFromCallable(gv,zlambda);
        Dune::PDELab::interpolate(zfunc,*disc.gfs,*coarse_vectors[3]);
      }

    for (int k = dim+1; k<nvec; k++)
    {
      for (size_t i = 0; i<istl_coarse_vectors[k]->N(); i++)
      {
        (*istl_coarse_vectors[k])[i] = (double)(rand() % 100);  
      }
    }


    // debug vectors
    VEC debug1(*disc.gfs);
    debug1 = 0.0;
    ISTLV& istldebug1 = Dune::PDELab::Backend::native(debug1);
    VEC debug2(*disc.gfs);
    debug2 = 0.0;
    ISTLV& istldebug2 = Dune::PDELab::Backend::native(debug2);
    VEC debug3(*disc.gfs);
    debug3 = 0.0;
    ISTLV& istldebug3 = Dune::PDELab::Backend::native(debug3);
    
    // set up an ISTL solver and solve
    bool twolevel = ptree.get("geneo.twolevel",(bool)true);
    Dune::MatrixAdapter<ISTLM,ISTLV,ISTLV> op(istlA);
    //SimpleAdditiveSchwarz<ISTLM,ISTLV> prec(matrices,flsdi.local_to_global,boundary_dofs);
    //RestrictedAdditiveSchwarz<ISTLM,ISTLV> prec(matrices,flsdi.local_to_global,flsdi.global_to_local,boundary_dofs,pum);
    std::cout << "constructing preconditioner" << std::endl;
    // TwoLevelAdditiveSchwarz<ISTLM,ISTLV> prec(pistlA,istl_coarse_vectors,pistldirichlet,floating,matrices, flsdi.local_to_global,flsdi.global_to_local,boundary_dofs,pum,twolevel, debug1,debug2,debug3);
    TwoLevelGenEO<ISTLM,ISTLV> prec(pistlA,istl_coarse_vectors,pistldirichlet,floating,matrices, flsdi.local_to_global,flsdi.global_to_local,boundary_dofs,pum,twolevel, debug1,debug2,debug3,ptree);
    //Dune::SeqILU<ISTLM,ISTLV,ISTLV> prec(istlA,1.0);
    //Dune::Richardson<ISTLV,ISTLV> prec;
    Dune::CGSolver<ISTLV> solver(op,prec,1e-8,ptree.get<int>("solver.maximum_iterations",1000),2,1);
    // Dune::RestartedGMResSolver<ISTLV> solver(op,prec,1e-8,70,ptree.get<int>("solver.maximum_iterations",1000),2);
    Dune::InverseOperatorResult stat;
    ISTLV istlv(istlx);
    istlv = (double)0.0;
    solver.apply(istlv,istlr,stat);
    istlx -= istlv;


    // std::vector<VEC> local_coarse_basis_output(10,debug1);
    // std::vector<std::vector<VEC>> coarse_basis_output(subdomains,local_coarse_basis_output);

    // // write some coarse basis vectors into debug output
    // for (int p = 0; p<subdomains; p++)
    // {
    //   for (int k = 0; k<ptree.get<int>("geneo.number_of_eigenvectors",1); k++)
    //   {
    //     ISTLV& temporary_istl = Dune::PDELab::Backend::native(coarse_basis_output[p][k]);
    //     temporary_istl = prec.prolongated_basis_vector(p,k);
    //   }
    // }

    // istldebug1 = prec.prolongated_basis_vector(7,7);   
    // istldebug2 = prec.prolongated_basis_vector(7,8);
    // istldebug3 = prec.prolongated_basis_vector(7,9);

    //===================
    // visualization
    //==================
    
    // visualize result
    typename Discretization::DGF dgf(*disc.gfs,x);
    typename Discretization::DGF dgfpu(*disc.gfs,puvec);
    typename Discretization::DGF dgfdebug1(*disc.gfs,debug1);
    typename Discretization::DGF dgfdebug2(*disc.gfs,debug2);
    typename Discretization::DGF dgfdebug3(*disc.gfs,debug3);
    std::string outfile = ptree.get("output.filename","output");
    Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::nonconforming);
    vtkwriter.addCellData(partition,"partition");
    vtkwriter.addCellData(adomain,"adomain");
    vtkwriter.addCellData(k0,"k0");
    using VTKF = Dune::PDELab::VTKGridFunctionAdapter<decltype(dgf)>;
    vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(dgf,"fesol")));
    vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(dgfpu,"pu")));
    vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(dgfdebug1,"debug1")));
    vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(dgfdebug2,"debug2")));
    vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(dgfdebug3,"debug3")));
    // for (int p = 0; p<subdomains; p++)
    // {
    //   for (int k = 0; k<ptree.get<int>("geneo.number_of_eigenvectors",10); k++)
    //   {
    //     auto dgf_coarse_basis_output = std::make_shared<Discretization::DGF>(*disc.gfs,coarse_basis_output[p][k]);
    //     vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(dgf_coarse_basis_output,"coarse_basis["+std::to_string(p)+"]["+std::to_string(k)+"]")));
    //   }
    // }

    vtkwriter.write(outfile,Dune::VTK::appendedraw);
#else
    std::cout << "Needs UGGrid" << std::endl;
#endif
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
