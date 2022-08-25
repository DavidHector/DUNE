#ifndef DUNE_MULTIDOMAINGRID_HH
#define DUNE_MULTIDOMAINGRID_HH

#include <dune/common/parallel/communication.hh>
#include <dune/grid/multidomaingrid/multidomaingrid.hh>
#include <dune/grid/multidomaingrid/multidomainmcmgmapper.hh>

namespace Dune {

using mdgrid::MultiDomainGrid;
using mdgrid::MultiDomainMCMGMapper;


namespace Capabilities {

  template<typename Grid>
  struct isAdaptable
  {
    static const bool v = true;
  };

}


// MultiDomainGrid capabilities

namespace Capabilities {

  template<class HostGrid, typename MDGridTraits, int codim>
  struct hasEntity<MultiDomainGrid<HostGrid,MDGridTraits>, codim>
  {
    static const bool v = hasEntity<HostGrid,codim>::v;
  };

  template<class HostGrid, typename MDGridTraits, int codim>
  struct hasEntityIterator<MultiDomainGrid<HostGrid,MDGridTraits>, codim>
  {
    static const bool v = hasEntityIterator<HostGrid,codim>::v;
  };

  template<class HostGrid, typename MDGridTraits, int codim>
  struct canCommunicate<MultiDomainGrid<HostGrid,MDGridTraits>, codim>
  {
    static const bool v = canCommunicate<HostGrid,codim>::v;
  };


  template<class HostGrid, typename MDGridTraits>
  struct isLevelwiseConforming<MultiDomainGrid<HostGrid,MDGridTraits> >
  {
    static const bool v = isLevelwiseConforming<HostGrid>::v;
  };


  template<class HostGrid, typename MDGridTraits>
  struct isLeafwiseConforming<MultiDomainGrid<HostGrid,MDGridTraits> >
  {
    static const bool v = isLeafwiseConforming<HostGrid>::v;
  };


  template<class HostGrid, typename MDGridTraits>
  struct hasBackupRestoreFacilities<MultiDomainGrid<HostGrid,MDGridTraits> >
  {
    static const bool v = false;
  };


  template<typename HostGrid, typename MDGridTraits>
  struct threadSafe<MultiDomainGrid<HostGrid,MDGridTraits> >
  {
    static const bool v = false;
  };


  template<typename HostGrid, typename MDGridTraits>
  struct viewThreadSafe<MultiDomainGrid<HostGrid,MDGridTraits> >
  {
    static const bool v = viewThreadSafe<HostGrid>::v;
  };


  template<typename HostGrid, typename MDGridTraits>
  struct hasSingleGeometryType<MultiDomainGrid<HostGrid,MDGridTraits> >
  {
    static const bool v = hasSingleGeometryType<HostGrid>::v;
    static const unsigned int topologyId = hasSingleGeometryType<HostGrid>::topologyId;
  };


  template<typename HostGrid, typename MDGridTraits, int codim>
  struct hasGeometry<MultiDomainGrid<HostGrid,MDGridTraits>,codim>
  {
    static const bool v = hasGeometry<HostGrid,codim>::v;
  };

} // namespace Capabilities



// SubDomainGrid capabilities

namespace Capabilities {

  template<class MDGrid, int codim>
  struct hasEntity<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid>, codim>
  {
    static const bool v = hasEntity<MDGrid,codim>::v;
  };

  template<class MDGrid, int codim>
  struct hasEntityIterator<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid>, codim>
  {
    static const bool v = hasEntityIterator<MDGrid,codim>::v;
  };

  template<class MDGrid, int codim>
  struct canCommunicate<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid>, codim>
  {
    static const bool v = canCommunicate<MDGrid,codim>::v;
  };


  template<class MDGrid>
  struct isLevelwiseConforming<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid> >
  {
    static const bool v = isLevelwiseConforming<MDGrid>::v;
  };


  template<class MDGrid>
  struct isLeafwiseConforming<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid> >
  {
    static const bool v = isLeafwiseConforming<MDGrid>::v;
  };


  template<class MDGrid>
  struct hasBackupRestoreFacilities<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid> >
  {
    static const bool v = false;
  };


  template<typename MDGrid>
  struct threadSafe<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid> >
  {
    static const bool v = false;
  };


  template<typename MDGrid>
  struct viewThreadSafe<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid> >
  {
    static const bool v = viewThreadSafe<MDGrid>::v;
  };


  template<typename MDGrid>
  struct hasSingleGeometryType<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid> >
  {
    static const bool v = hasSingleGeometryType<MDGrid>::v;
    static const unsigned int topologyId = hasSingleGeometryType<MDGrid>::topologyId;
  };

  template<typename MDGrid, int codim>
  struct hasGeometry<Dune::mdgrid::subdomain::SubDomainGrid<MDGrid>,codim >
  {
    static const bool v = hasGeometry<MDGrid,codim>::v;
  };


} // namespace Capabilities

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_HH
