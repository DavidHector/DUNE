// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_HYDRO_NONLINEARDIFFUSION_VELOCITY_HH
#define DUNE_HYDRO_NONLINEARDIFFUSION_VELOCITY_HH

// dune-core includes
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/typetraits.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/localfunctions/raviartthomas/raviartthomascube.hh>

// dune-pdelab includes
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/convectiondiffusionparameter.hh>


/** \brief Provides a vector-valued grid function for velocity field

    using RT0 interpolation on a rectangular cell

    - Model : implements parameter class
    - GFS   : P0 grid function space on given mesh
    - Z     : P0 coefficient vector
*/
template<typename  Model, typename GFS, typename Z>
class NonlinearDiffusionFVVelocity
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename GFS::Traits::GridView,
             typename GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
             GFS::Traits::GridView::dimension,
             Dune::FieldVector<typename GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,GFS::Traits::GridView::dimension> >,
             NonlinearDiffusionFVVelocity<Model,GFS,Z> >
{
  // extract useful types
  typedef typename GFS::Traits::GridView GV;
  enum { dim = GV::dimension };
  typedef typename GV::IndexSet IndexSet;
  typedef typename GV::Grid::ctype DF;
  typedef typename GFS::Traits::FiniteElementType FiniteElementType;
  typedef typename FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
  typedef typename FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;
  typedef typename Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0>::Traits::LocalBasisType::Traits::RangeType RT0RangeType;

  // member variables
  const Model& model;
  GFS gfs;
  GV gv;
  const IndexSet& is;
  const Z& solution;
  const Z& bathymmetry;

  Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0> rt0fe;
  mutable std::vector<RT0RangeType> rt0vectors;
  mutable Dune::FieldMatrix<DF,dim,dim> B;
  mutable RF determinant;
  mutable int cachedindex;
  RF time;
  RF residualheight;

  typedef Dune::FieldVector<RF,2*dim> RT0Coeffs;
  std::vector<RT0Coeffs> storedcoeffs;

  // A DataHandle class to exchange RT0 coefficients in the overlap
  template<class GV, class V> // mapper type and vector type
  class VectorExchange
    : public Dune::CommDataHandleIF<VectorExchange<GV,V>,
                                    typename V::value_type>
  {
    typedef typename GV::IndexSet IndexSet;
    typedef typename IndexSet::IndexType IndexType;

    GV gv;
    V& c;
    const IndexSet& indexSet;

  public:
    //! export type of data for message buffer
    typedef typename V::value_type DataType;

    //! constructor
    VectorExchange (const GV& gv_, V& c_)
      : gv(gv_), c(c_), indexSet(gv.indexSet())
    {}

    //! returns true if data for this codim should be communicated
    bool contains (int dim, int codim) const
    {
      return (codim==0);
    }

    //! returns true if size per entity of given dim and codim is a constant
    bool fixedSize (int dim, int codim) const
    {
      return true;
    }

    /*! how many objects of type DataType have to be sent for a given entity

      Note: Only the sender side needs to know this size.
    */
    template<class EntityType>
    size_t size (const EntityType& e) const
    {
      return 1;
    }

    //! pack data from user to message buffer
    template<class MessageBuffer, class EntityType>
    void gather (MessageBuffer& buff, const EntityType& e) const
    {
      buff.write(c[indexSet.index(e)]);
    }

    /*! unpack data from message buffer to user

      n is the number of objects sent by the sender
    */
    template<class MessageBuffer, class EntityType>
    void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
    {
      DataType x;
      buff.read(x);
      c[indexSet.index(e)]=x;
    }
  };


public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Traits,NonlinearDiffusionFVVelocity<Model,GFS,Z> > BaseT;

  NonlinearDiffusionFVVelocity (const Model& model_, const GFS& gfs_, const Z& solution_, const Z& bathymmetry_, RF residualheight_=0.0)
    : model(model_), gfs(gfs_), gv(gfs.gridView()), is(gv.indexSet()), solution(solution_), bathymmetry(bathymmetry_),
      rt0fe(), rt0vectors(rt0fe.localBasis().size()), cachedindex(-1), time(0.0), storedcoeffs(is.size(0)), residualheight(residualheight_)
  {
    update();
  }

  void update ()
  {
    // update invalidates cached data
    cachedindex = -1;
    
    // compute RT0 coefficients for all interior cells
    for (const auto& cell_inside : elements(gv,Dune::Partitions::interior))
      {
        // evaluate inside cell
        int index_inside = is.index(cell_inside);
        auto insidegeo = cell_inside.geometry();
        auto inside_global = insidegeo.center();
        auto inside_cell_center_local = referenceElement(insidegeo).position(0,0);

        // evaluate cell data
        auto u_inside = Dune::PDELab::Backend::native(solution)[index_inside][0];
        auto b_inside = Dune::PDELab::Backend::native(bathymmetry)[index_inside][0];
        auto k_inside = model.k(cell_inside,inside_cell_center_local);

        // for coefficient computation
        RF vn[2*dim];    // normal velocities
        B = cell_inside.geometry().jacobianInverseTransposed(inside_cell_center_local); // the transformation. Assume it is linear
        determinant = B.determinant();

        // loop over cell neighbors
        for (const auto& intersection : intersections(gv,cell_inside))
          {
            // set to zero for processor boundary
            vn[intersection.indexInInside()] = 0.0;

            // face geometry
            auto facegeo = intersection.geometry();

            // interior face
            if (intersection.neighbor())
              {
                // evaluate outside cell
                const auto& cell_outside = intersection.outside();
                int index_outside = is.index(cell_outside);
                auto outsidegeo = cell_outside.geometry();
                auto outside_global = outsidegeo.center();
                const auto outside_cell_center_local = referenceElement(outsidegeo).position(0,0);

                // distance between the two cell centers
                outside_global -= inside_global;
                auto distance = outside_global.two_norm();

                // evaluate data
                auto u_outside = Dune::PDELab::Backend::native(solution)[index_outside][0];
                auto b_outside = Dune::PDELab::Backend::native(bathymmetry)[index_outside][0];
                auto k_outside = model.k(cell_outside,outside_cell_center_local);

                if (u_inside-b_inside<=residualheight ||  u_outside-b_outside<=residualheight)
                  vn[intersection.indexInInside()] = 0.0;
                else
                  {
                    // gradient computed from current solution
                    auto VN = (u_inside-u_outside)/distance;

                    // permeability
                    auto k = 2.0*(k_inside+1e-15)*(k_outside+1e-15)/(k_inside+k_outside+2e-15);
    
                    // nonlinearity
                    auto phi = model.phi(b_inside,b_outside,u_inside,u_outside,VN);
    
                    // set coefficient
                    vn[intersection.indexInInside()] = k*phi*VN;
                  }
              }

            // boundary face
            if (intersection.boundary())
              {
                // check for Dirichlet boundary condition
                auto facecenterlocal = referenceElement(facegeo).position(0,0);
                bool isdirichlet = model.b(intersection,facecenterlocal);

                if (isdirichlet)
                  {
                    // Dirichlet boundary condition
                    auto outside_global = facegeo.center();
                    outside_global -= inside_global;
                    auto distance = outside_global.two_norm();

                    auto facecenterinelement=intersection.geometryInInside().center();
                    auto u_outside = model.g(cell_inside,facecenterinelement);
                    
                    // gradient computed from current solution
                    auto VN = (u_inside-u_outside)/distance;

                    vn[intersection.indexInInside()] = k_inside*VN;
                  }
                else
                  {
                    // Neumann boundary condition
                    vn[intersection.indexInInside()] = model.j(intersection,facecenterlocal);
                  }
              }

            // compute coefficient
            Dune::FieldVector<DF,dim> vstar=intersection.centerUnitOuterNormal(); // normal on tranformef element
            vstar *= vn[intersection.indexInInside()];
            Dune::FieldVector<RF,dim> normalhat(0); // normal on reference element
            if (intersection.indexInInside()%2==0)
              normalhat[intersection.indexInInside()/2] = -1.0;
            else
              normalhat[intersection.indexInInside()/2] =  1.0;
            Dune::FieldVector<DF,dim> vstarhat(0);
            B.umtv(vstar,vstarhat); // Piola backward transformation
            vstarhat *= determinant;
            storedcoeffs[index_inside][intersection.indexInInside()] = vstarhat*normalhat;
          }
      }

    // communicate coefficients in overlap
    VectorExchange<GV,std::vector<RT0Coeffs> > dh(gv,storedcoeffs);
    if (gv.grid().comm().size()>1)
      gv.grid().communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  }

  // set time where operator is to be evaluated (i.e. end of the time intervall)
  void set_time (RF time_)
  {
    time = time_;
  }

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    // local cell number
    int index = is.index(e);

    // compute velocity on reference element
    rt0fe.localBasis().evaluateFunction(x,rt0vectors);
    typename Traits::RangeType yhat(0);
    for (unsigned int i=0; i<rt0fe.localBasis().size(); i++)
      yhat.axpy(storedcoeffs[index][i],rt0vectors[i]);

    // apply Piola transformation
    if (index != cachedindex)
      {
        B = e.geometry().jacobianTransposed(x); // the transformation. Assume it is linear
        determinant = B.determinant();
        cachedindex = index;
      }
    y = 0;
    B.umtv(yhat,y);
    y *= determinant;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }
};


/** \brief Provides a vector-valued grid function for velocity field

    using RT0 interpolation on a rectangular cell

    - Model : implements parameter class
    - GFS   : P0 grid function space on given mesh
    - Z     : P0 coefficient vector
*/
template<typename  Model, typename GFS, typename Z>
class NonlinearDiffusionFVVelocitySurface
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename GFS::Traits::GridView,
             typename GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
             GFS::Traits::GridView::dimension,
             Dune::FieldVector<typename GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,GFS::Traits::GridView::dimension> >,
             NonlinearDiffusionFVVelocitySurface<Model,GFS,Z> >
{
  // extract useful types
  typedef typename GFS::Traits::GridView GV;
  enum { dim = GV::dimension };
  typedef typename GV::IndexSet IndexSet;
  typedef typename GV::Grid::ctype DF;
  typedef typename GFS::Traits::FiniteElementType FiniteElementType;
  typedef typename FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
  typedef typename FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;
  typedef typename Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0>::Traits::LocalBasisType::Traits::RangeType RT0RangeType;

  // member variables
  const Model& model;
  GFS gfs;
  GV gv;
  const IndexSet& is;
  const Z& solution;
  const Z& bathymmetry;

  Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0> rt0fe;
  mutable std::vector<RT0RangeType> rt0vectors;
  mutable Dune::FieldMatrix<DF,dim,dim> B;
  mutable RF determinant;
  mutable int cachedindex;
  RF time;
  RF residualheight;

  typedef Dune::FieldVector<RF,2*dim> RT0Coeffs;
  std::vector<RT0Coeffs> storedcoeffs;

  // A DataHandle class to exchange RT0 coefficients in the overlap
  template<class GV, class V> // mapper type and vector type
  class VectorExchange
    : public Dune::CommDataHandleIF<VectorExchange<GV,V>,
                                    typename V::value_type>
  {
    typedef typename GV::IndexSet IndexSet;
    typedef typename IndexSet::IndexType IndexType;

    GV gv;
    V& c;
    const IndexSet& indexSet;

  public:
    //! export type of data for message buffer
    typedef typename V::value_type DataType;

    //! constructor
    VectorExchange (const GV& gv_, V& c_)
      : gv(gv_), c(c_), indexSet(gv.indexSet())
    {}

    //! returns true if data for this codim should be communicated
    bool contains (int dim, int codim) const
    {
      return (codim==0);
    }

    //! returns true if size per entity of given dim and codim is a constant
    bool fixedSize (int dim, int codim) const
    {
      return true;
    }

    /*! how many objects of type DataType have to be sent for a given entity

      Note: Only the sender side needs to know this size.
    */
    template<class EntityType>
    size_t size (const EntityType& e) const
    {
      return 1;
    }

    //! pack data from user to message buffer
    template<class MessageBuffer, class EntityType>
    void gather (MessageBuffer& buff, const EntityType& e) const
    {
      buff.write(c[indexSet.index(e)]);
    }

    /*! unpack data from message buffer to user

      n is the number of objects sent by the sender
    */
    template<class MessageBuffer, class EntityType>
    void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
    {
      DataType x;
      buff.read(x);
      c[indexSet.index(e)]=x;
    }
  };


public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Traits,NonlinearDiffusionFVVelocitySurface<Model,GFS,Z> > BaseT;

  NonlinearDiffusionFVVelocitySurface (const Model& model_, const GFS& gfs_, const Z& solution_, const Z& bathymmetry_, RF residualheight_=0.0)
    : model(model_), gfs(gfs_), gv(gfs.gridView()), is(gv.indexSet()), solution(solution_), bathymmetry(bathymmetry_),
      rt0fe(), rt0vectors(rt0fe.localBasis().size()), cachedindex(-1), time(0.0), storedcoeffs(is.size(0)), residualheight(residualheight_)
  {
    update();
  }

  void update ()
  {
    // update invalidates cached data
    cachedindex = -1;
    
    // compute RT0 coefficients for all interior cells
    for (const auto& cell_inside : elements(gv,Dune::Partitions::interior))
      {
        // evaluate inside cell
        int index_inside = is.index(cell_inside);
        auto insidegeo = cell_inside.geometry();
        auto inside_global = insidegeo.center();
        auto inside_cell_center_local = referenceElement(insidegeo).position(0,0);

        // evaluate cell data
        auto u_inside = Dune::PDELab::Backend::native(solution)[index_inside][0];
        auto b_inside = Dune::PDELab::Backend::native(bathymmetry)[index_inside][0];
        auto k_inside = model.k_surface(cell_inside,inside_cell_center_local);

        // for coefficient computation
        RF vn[2*dim];    // normal velocities
        B = cell_inside.geometry().jacobianInverseTransposed(inside_cell_center_local); // the transformation. Assume it is linear
        determinant = B.determinant();

        // loop over cell neighbors
        for (const auto& intersection : intersections(gv,cell_inside))
          {
            // set to zero for processor boundary
            vn[intersection.indexInInside()] = 0.0;

            // face geometry
            auto facegeo = intersection.geometry();

            // interior face
            if (intersection.neighbor())
              {
                // evaluate outside cell
                const auto& cell_outside = intersection.outside();
                int index_outside = is.index(cell_outside);
                auto outsidegeo = cell_outside.geometry();
                auto outside_global = outsidegeo.center();
                const auto outside_cell_center_local = referenceElement(outsidegeo).position(0,0);

                // distance between the two cell centers
                outside_global -= inside_global;
                auto distance = outside_global.two_norm();

                // evaluate data
                auto u_outside = Dune::PDELab::Backend::native(solution)[index_outside][0];
                auto b_outside = Dune::PDELab::Backend::native(bathymmetry)[index_outside][0];
                auto k_outside = model.k_surface(cell_outside,outside_cell_center_local);

                if (u_inside-b_inside<=residualheight ||  u_outside-b_outside<=residualheight)
                  vn[intersection.indexInInside()] = 0.0;
                else
                  {
                    // gradient computed from current solution
                    auto VN = (u_inside-u_outside)/distance;

                    // permeability
                    auto k = 2.0*(k_inside+1e-15)*(k_outside+1e-15)/(k_inside+k_outside+2e-15);
    
                    // nonlinearity
                    auto phi = model.phi_surface(b_inside,b_outside,u_inside,u_outside,VN);
    
                    // set coefficient
                    vn[intersection.indexInInside()] = k*phi*VN;
                  }
              }

            // boundary face
            if (intersection.boundary())
              {
                // check for Dirichlet boundary condition
                auto facecenterlocal = referenceElement(facegeo).position(0,0);
                bool isdirichlet = model.b_surface(intersection,facecenterlocal);

                if (isdirichlet)
                  {
                    // Dirichlet boundary condition
                    auto outside_global = facegeo.center();
                    outside_global -= inside_global;
                    auto distance = outside_global.two_norm();

                    auto facecenterinelement=intersection.geometryInInside().center();
                    auto u_outside = model.g_surface(cell_inside,facecenterinelement);
                    
                    // gradient computed from current solution
                    auto VN = (u_inside-u_outside)/distance;

                    vn[intersection.indexInInside()] = k_inside*VN;
                  }
                else
                  {
                    // Neumann boundary condition
                    vn[intersection.indexInInside()] = model.j_surface(intersection,facecenterlocal);
                  }
              }

            // compute coefficient
            Dune::FieldVector<DF,dim> vstar=intersection.centerUnitOuterNormal(); // normal on tranformef element
            vstar *= vn[intersection.indexInInside()];
            Dune::FieldVector<RF,dim> normalhat(0); // normal on reference element
            if (intersection.indexInInside()%2==0)
              normalhat[intersection.indexInInside()/2] = -1.0;
            else
              normalhat[intersection.indexInInside()/2] =  1.0;
            Dune::FieldVector<DF,dim> vstarhat(0);
            B.umtv(vstar,vstarhat); // Piola backward transformation
            vstarhat *= determinant;
            storedcoeffs[index_inside][intersection.indexInInside()] = vstarhat*normalhat;
          }
      }

    // communicate coefficients in overlap
    VectorExchange<GV,std::vector<RT0Coeffs> > dh(gv,storedcoeffs);
    if (gv.grid().comm().size()>1)
      gv.grid().communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  }

  // set time where operator is to be evaluated (i.e. end of the time intervall)
  void set_time (RF time_)
  {
    time = time_;
  }

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    // local cell number
    int index = is.index(e);

    // compute velocity on reference element
    rt0fe.localBasis().evaluateFunction(x,rt0vectors);
    typename Traits::RangeType yhat(0);
    for (unsigned int i=0; i<rt0fe.localBasis().size(); i++)
      yhat.axpy(storedcoeffs[index][i],rt0vectors[i]);

    // apply Piola transformation
    if (index != cachedindex)
      {
        B = e.geometry().jacobianTransposed(x); // the transformation. Assume it is linear
        determinant = B.determinant();
        cachedindex = index;
      }
    y = 0;
    B.umtv(yhat,y);
    y *= determinant;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }
};



/** \brief Provides a vector-valued grid function for velocity field

    using RT0 interpolation on a rectangular cell

    - Model : implements parameter class
    - GFS   : P0 grid function space on given mesh
    - Z     : P0 coefficient vector
*/
template<typename  Model, typename GFS, typename Z>
class NonlinearDiffusionFVVelocityGroundwater
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename GFS::Traits::GridView,
             typename GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,
             GFS::Traits::GridView::dimension,
             Dune::FieldVector<typename GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType,GFS::Traits::GridView::dimension> >,
             NonlinearDiffusionFVVelocityGroundwater<Model,GFS,Z> >
{
  // extract useful types
  typedef typename GFS::Traits::GridView GV;
  enum { dim = GV::dimension };
  typedef typename GV::IndexSet IndexSet;
  typedef typename GV::Grid::ctype DF;
  typedef typename GFS::Traits::FiniteElementType FiniteElementType;
  typedef typename FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
  typedef typename FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;
  typedef typename Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0>::Traits::LocalBasisType::Traits::RangeType RT0RangeType;

  // member variables
  const Model& model;
  GFS gfs;
  GV gv;
  const IndexSet& is;
  const Z& solution;
  const Z& bathymmetry;

  Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0> rt0fe;
  mutable std::vector<RT0RangeType> rt0vectors;
  mutable Dune::FieldMatrix<DF,dim,dim> B;
  mutable RF determinant;
  mutable int cachedindex;
  RF time;
  RF residualheight;

  typedef Dune::FieldVector<RF,2*dim> RT0Coeffs;
  std::vector<RT0Coeffs> storedcoeffs;

  // A DataHandle class to exchange RT0 coefficients in the overlap
  template<class GV, class V> // mapper type and vector type
  class VectorExchange
    : public Dune::CommDataHandleIF<VectorExchange<GV,V>,
                                    typename V::value_type>
  {
    typedef typename GV::IndexSet IndexSet;
    typedef typename IndexSet::IndexType IndexType;

    GV gv;
    V& c;
    const IndexSet& indexSet;

  public:
    //! export type of data for message buffer
    typedef typename V::value_type DataType;

    //! constructor
    VectorExchange (const GV& gv_, V& c_)
      : gv(gv_), c(c_), indexSet(gv.indexSet())
    {}

    //! returns true if data for this codim should be communicated
    bool contains (int dim, int codim) const
    {
      return (codim==0);
    }

    //! returns true if size per entity of given dim and codim is a constant
    bool fixedSize (int dim, int codim) const
    {
      return true;
    }

    /*! how many objects of type DataType have to be sent for a given entity

      Note: Only the sender side needs to know this size.
    */
    template<class EntityType>
    size_t size (const EntityType& e) const
    {
      return 1;
    }

    //! pack data from user to message buffer
    template<class MessageBuffer, class EntityType>
    void gather (MessageBuffer& buff, const EntityType& e) const
    {
      buff.write(c[indexSet.index(e)]);
    }

    /*! unpack data from message buffer to user

      n is the number of objects sent by the sender
    */
    template<class MessageBuffer, class EntityType>
    void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
    {
      DataType x;
      buff.read(x);
      c[indexSet.index(e)]=x;
    }
  };


public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> > Traits;
  typedef Dune::PDELab::GridFunctionBase<Traits,NonlinearDiffusionFVVelocityGroundwater<Model,GFS,Z> > BaseT;

  NonlinearDiffusionFVVelocityGroundwater (const Model& model_, const GFS& gfs_, const Z& solution_, const Z& bathymmetry_, RF residualheight_=0.0)
    : model(model_), gfs(gfs_), gv(gfs.gridView()), is(gv.indexSet()), solution(solution_), bathymmetry(bathymmetry_),
      rt0fe(), rt0vectors(rt0fe.localBasis().size()), cachedindex(-1), time(0.0), storedcoeffs(is.size(0)), residualheight(residualheight_)
  {
    update();
  }

  void update ()
  {
    // update invalidates cached data
    cachedindex = -1;
    
    // compute RT0 coefficients for all interior cells
    for (const auto& cell_inside : elements(gv,Dune::Partitions::interior))
      {
        // evaluate inside cell
        int index_inside = is.index(cell_inside);
        auto insidegeo = cell_inside.geometry();
        auto inside_global = insidegeo.center();
        auto inside_cell_center_local = referenceElement(insidegeo).position(0,0);

        // evaluate cell data
        auto u_inside = Dune::PDELab::Backend::native(solution)[index_inside][0];
        auto b_inside = Dune::PDELab::Backend::native(bathymmetry)[index_inside][0];
        auto k_inside = model.k_groundwater(cell_inside,inside_cell_center_local);

        // for coefficient computation
        RF vn[2*dim];    // normal velocities
        B = cell_inside.geometry().jacobianInverseTransposed(inside_cell_center_local); // the transformation. Assume it is linear
        determinant = B.determinant();

        // loop over cell neighbors
        for (const auto& intersection : intersections(gv,cell_inside))
          {
            // set to zero for processor boundary
            vn[intersection.indexInInside()] = 0.0;

            // face geometry
            auto facegeo = intersection.geometry();

            // interior face
            if (intersection.neighbor())
              {
                // evaluate outside cell
                const auto& cell_outside = intersection.outside();
                int index_outside = is.index(cell_outside);
                auto outsidegeo = cell_outside.geometry();
                auto outside_global = outsidegeo.center();
                const auto outside_cell_center_local = referenceElement(outsidegeo).position(0,0);

                // distance between the two cell centers
                outside_global -= inside_global;
                auto distance = outside_global.two_norm();

                // evaluate data
                auto u_outside = Dune::PDELab::Backend::native(solution)[index_outside][0];
                auto b_outside = Dune::PDELab::Backend::native(bathymmetry)[index_outside][0];
                auto k_outside = model.k_groundwater(cell_outside,outside_cell_center_local);

                if (u_inside-b_inside<=residualheight ||  u_outside-b_outside<=residualheight)
                  vn[intersection.indexInInside()] = 0.0;
                else
                  {
                    // gradient computed from current solution
                    auto VN = (u_inside-u_outside)/distance;

                    // permeability
                    auto k = 2.0*(k_inside+1e-15)*(k_outside+1e-15)/(k_inside+k_outside+2e-15);
    
                    // nonlinearity
                    auto phi = model.phi_groundwater(b_inside,b_outside,u_inside,u_outside,VN);
    
                    // set coefficient
                    vn[intersection.indexInInside()] = k*phi*VN;
                  }
              }

            // boundary face
            if (intersection.boundary())
              {
                // check for Dirichlet boundary condition
                auto facecenterlocal = referenceElement(facegeo).position(0,0);
                bool isdirichlet = model.b_groundwater(intersection,facecenterlocal);

                if (isdirichlet)
                  {
                    // Dirichlet boundary condition
                    auto outside_global = facegeo.center();
                    outside_global -= inside_global;
                    auto distance = outside_global.two_norm();

                    auto facecenterinelement=intersection.geometryInInside().center();
                    auto u_outside = model.g_groundwater(cell_inside,facecenterinelement);
                    
                    // gradient computed from current solution
                    auto VN = (u_inside-u_outside)/distance;

                    vn[intersection.indexInInside()] = k_inside*VN;
                  }
                else
                  {
                    // Neumann boundary condition
                    vn[intersection.indexInInside()] = model.j_groundwater(intersection,facecenterlocal);
                  }
              }

            // compute coefficient
            Dune::FieldVector<DF,dim> vstar=intersection.centerUnitOuterNormal(); // normal on tranformef element
            vstar *= vn[intersection.indexInInside()];
            Dune::FieldVector<RF,dim> normalhat(0); // normal on reference element
            if (intersection.indexInInside()%2==0)
              normalhat[intersection.indexInInside()/2] = -1.0;
            else
              normalhat[intersection.indexInInside()/2] =  1.0;
            Dune::FieldVector<DF,dim> vstarhat(0);
            B.umtv(vstar,vstarhat); // Piola backward transformation
            vstarhat *= determinant;
            storedcoeffs[index_inside][intersection.indexInInside()] = vstarhat*normalhat;
          }
      }

    // communicate coefficients in overlap
    VectorExchange<GV,std::vector<RT0Coeffs> > dh(gv,storedcoeffs);
    if (gv.grid().comm().size()>1)
      gv.grid().communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  }

  // set time where operator is to be evaluated (i.e. end of the time intervall)
  void set_time (RF time_)
  {
    time = time_;
  }

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    // local cell number
    int index = is.index(e);

    // compute velocity on reference element
    rt0fe.localBasis().evaluateFunction(x,rt0vectors);
    typename Traits::RangeType yhat(0);
    for (unsigned int i=0; i<rt0fe.localBasis().size(); i++)
      yhat.axpy(storedcoeffs[index][i],rt0vectors[i]);

    // apply Piola transformation
    if (index != cachedindex)
      {
        B = e.geometry().jacobianTransposed(x); // the transformation. Assume it is linear
        determinant = B.determinant();
        cachedindex = index;
      }
    y = 0;
    B.umtv(yhat,y);
    y *= determinant;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }
};



/** \brief compute maximum velocity from vector-valued grid function
 */
template<typename T>
typename T::Traits::RangeFieldType maxvelocity (const T& t)
{
  // extract useful types
  typedef typename T::Traits::GridViewType GV;
  typedef typename T::Traits::RangeFieldType RF;
  typedef typename T::Traits::DomainType DomainType;
  typedef typename T::Traits::RangeType RangeType;

  // variables
  const int dim = T::Traits::GridViewType::dimension;
  RF maximum = 0.0;
  GV gv(t.getGridView());

  // loop over the grid
  for (const auto& cell : elements(gv,Dune::Partitions::interior))
    {
      for (int d=0; d<dim; d++)
        {
          DomainType x0(0.5);
          x0[d] = 0.0;
          RangeType y0;
          t.evaluate(cell,x0,y0);
          maximum = std::max(maximum,y0.two_norm());
          DomainType x1(0.5);
          x1[d] = 1.0;
          RangeType y1;
          t.evaluate(cell,x1,y1);
          maximum = std::max(maximum,y1.two_norm());
        }
    }
  maximum = gv.comm().max(maximum);
  return maximum;
}


#endif
