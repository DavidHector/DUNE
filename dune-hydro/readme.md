# dune-hydro Module

## Availability

The latest version of dune-hydro is available from

    git clone https://parcomp-git.iwr.uni-heidelberg.de/peter/dune-hydro.git

## Introduction

The dune-hydro module implements several important models for hydrospheric flow and
transport. Initially, its focus will be on shallow flow models for surface and subsurface flows,
i.e. models where the dependence on depth is not resolved.

The *Version 1.0* covers the following models and methods:

- shallow groundwater flow
- shallow surface flow using the diffusive wave approximation
- transport of a dissolved substance of these flows
- discretization with cell-centered finite volumes and two-point flux on structured grids

In particular, there is no coupled surface and subsurface flow model yet.

An accompanying description of the models can be dowloaded as a git repository from

    git clone https://parcomp-git.iwr.uni-heidelberg.de/peter/shallow-hydro-models.git

## Installation 

dune-hydro is based on [dune-pdelab](https://www.dune-project.org/modules/dune-pdelab/).

It requires the installation of the [GDAL Library](https://gdal.org/) (this dependence can currently not
be disabled) for reading [geo tiff files](https://en.wikipedia.org/wiki/GeoTIFF).
Such data can be be obtained for example via the [USGS earth explorer](https://earthexplorer.usgs.gov/).

## Code Structure

###src

Contains several examples:

- firstexample: A simple example for shallow groundwater flow with constant bathymmetry and a single sink centred in the middle pumping out water.
-  secondexample: A hump of water flowing down an inclined plane and making a pool.
-  thirdexample: flow in a complicated bathymmetry with constant rain. The region covers the Neckar valley from Heidelberg to Hirschhorn. 
-  fourthexample: the second example including concentration of a dissolved substance.

###dune/hydro

Contains discretisation schemes and supporting tools.


