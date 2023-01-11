#!/bin/bash
#
# This script can be used to download and install Dune with
# several external modules from source.
# Currently UG, HDF5 and FFTW3 are supported.
# If no argument is given it prints a list of software required for the build.

# here you may manually define command names
# comment out these definitions to overwrite the default from your system
#export F77=gfortran-mp-9
#export CC=gcc-mp-9
#export CXX=g++-mp-9
export CXXFLAGS='-I/opt/local/include -L/opt/local/lib -I/home/dhector/include -I/home/dhector/lib -O3 -DNDEBUG'
#export MPICC=mpicc-mpich-gcc9
#export MPICXX=mpicxx-mpich-gcc9
#export CMAKE=/opt/bin/cmake #im pool

# print requirements
print_requirements_message ()
{
    echo "Run this installer with "
    echo " "
    echo "   installer <directory>"
    echo " "
    echo "where <directory> is the name of a directory where to install to."
    echo "This directory is created if it does not exist."
    echo " "
    echo "In order to run this script you need the following software installed on your machine."
    echo "Use your package manager (or on the Mac macports or homebrew) to do this."
    echo " "
    echo "C++ Compiler. With pdelab gcc>=5.1 should be fine"
    echo "MPI-2 (including support for MPI IO)"
    echo "cmake>=3.0"
    echo "git"
    echo "wget"
    echo " "
    echo "Optionally you should consider to get:"
    echo "Paraview for visualization (http://www.paraview.org)"
    echo "Gmsh for mesh generation (http://gmsh.info)"
}

#
if [ $# -eq 0 ]
then
    print_requirements_message
    exit 1
fi

# usage: installer <directory>
#   directory: new directory in current directory where software is installed

# bash setup
set -x
set -e

# first we create a new directory, step into it and remember where we are
mkdir -p $1
cd $1
ROOT=$(pwd)

# make directory where to put external software
if [ ! "$EXTERNAL_HOME" ]; then
  EXTERNAL_HOME=$ROOT/external
  mkdir -p $EXTERNAL_HOME
fi

# extract default settings for the command names
if [ ! "$F77" ]; then
  F77=gfortran
fi
if [ ! "$CC" ]; then
CC=gcc
fi
if [ ! "$MPICC" ]; then
MPICC=mpicc
fi
if [ ! "$MPICXX" ]; then
MPICXX=mpicxx
fi
if [ ! "$CXX" ]; then
CXX=g++
fi
if [ ! "$CXXFLAGS" ]; then
CXXFLAGS="-O3 -DNDEBUG"
fi
CFLAGS="$CXXFLAGS"
if [ ! "$MAKE_FLAGS" ]; then
MAKE_FLAGS="-j2"
fi

# now lets define some functions

# install vector class library
install_vcl ()
{
    curl -O http://www.agner.org/optimize/vectorclass.zip
    unzip vectorclass.zip -d $EXTERNAL_HOME/vcl
}

# install arpack ng library
install_arpack-ng ()
{
    pushd $EXTERNAL_HOME
    git clone https://github.com/opencollab/arpack-ng.git arpack-src
    mkdir arpack-ng
    pushd arpack-src
    mkdir build
    pushd build
    cmake -DCMAKE_INSTALL_PREFIX=$EXTERNAL_HOME/arpack-ng -DCMAKE_C_COMPILER='gcc' -DCMAKE_CXX_COMPILER='g++' -DCMAKE_Fortran_COMPILER='gfortran' -D EXAMPLES=ON -D MPI=ON -D ICB=ON -D BUILD_SHARED_LIBS=ON ..
    make
    make install
    popd
    popd
    rm -rf arpack-src
    # now "install" arpackpp
    git clone https://github.com/m-reuter/arpackpp.git arpackpp-src
    pushd arpackpp-src
    mkdir build
    pushd build
    cmake -DCMAKE_INSTALL_PREFIX=$EXTERNAL_HOME/arpackpp -DCMAKE_C_COMPILER='gcc' -DCMAKE_CXX_COMPILER='g++' -DCMAKE_Fortran_COMPILER='gfortran' ..
    make install
    popd
    popd
    rm -rf arpackpp-src
    popd
}

# build fftw3
install_fftw3 ()
{
    FFTWVERSION=3.3.4
    pushd $EXTERNAL_HOME
    mkdir -p tarballs
    pushd tarballs
    wget http://www.fftw.org/fftw-$FFTWVERSION.tar.gz
    tar zxvf fftw-$FFTWVERSION.tar.gz
    pushd fftw-$FFTWVERSION
    ( ./configure CC=$CC MPICC=$MPICC CFLAGS="$CXXFLAGS" --enable-mpi --prefix=$EXTERNAL_HOME/fftw3 && make && make install) || exit $?
    popd
    rm -rf fftw-$FFTWVERSION
    popd
    popd
}

# build hdf5 with mpi support
install_hdf5 ()
{
    pushd $EXTERNAL_HOME
    mkdir -p tarballs
    pushd tarballs
#    HDF5VERSION=1.8.17
#    curl -O http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-$HDF5VERSION.tar.gz
    MAJOR=1.8
    MINOR=22
    HDF5VERSION=$MAJOR.$MINOR
    wget http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-$MAJOR/hdf5-$HDF5VERSION/src/hdf5-$HDF5VERSION.tar.gz
    tar zxvf hdf5-$HDF5VERSION.tar.gz
    pushd hdf5-$HDF5VERSION
    ( ./configure CC=$MPICC --enable-parallel --enable-static --prefix=$EXTERNAL_HOME/hdf5 && make && make install) || exit $?
    popd
#    rm -rf hdf5-$HDF5VERSION
    popd
    popd
}


# check out all the required dune modules
checkout_dune_core ()
{
    # dune core modules
    for i in common geometry grid istl localfunctions; do
	git clone https://gitlab.dune-project.org/core/dune-$i.git
	cd dune-$i
#	git checkout tags/v2.8.0rc1 -b 2.8.0rc1-branch
	git checkout releases/2.8
	cd ..
    done

    # modules in staging
    for i in functions uggrid typetree; do
	git clone https://gitlab.dune-project.org/staging/dune-$i.git
	cd dune-$i
#	git checkout tags/v2.8.0rc1 -b 2.8.0rc1-branch
	git checkout releases/2.8
	cd ..
    done

    # dune-alugrid
    git clone https://gitlab.dune-project.org/extensions/dune-alugrid.git
    cd dune-alugrid
    git checkout releases/2.8
    cd ..

    git clone https://gitlab.dune-project.org/extensions/dune-vtk.git
    cd dune-vtk 
    git checkout releases/2.8
    cd ..

	
    # install whitespace hook
    ./dune-common/bin/dunecontrol vcsetup
}

checkout_dune_pdelab ()
{
    # dune pdelab modules
	git clone https://gitlab.dune-project.org/pdelab/dune-pdelab.git
	cd dune-pdelab
#	git checkout releases/2.6
	cd ..

    # install whitespace hook
    ./dune-common/bin/dunecontrol vcsetup
}

checkout_dune_pdelab_applications ()
{
    # applications
#    git clone https://gitlab.dune-project.org/quality/dune-testtools.git
#    git clone https://gitlab.dune-project.org/extensions/dune-codegen.git
    git clone https://gitlab.dune-project.org/pdelab/dune-pdelab-tutorials.git
#    git clone https://gitlab.dune-project.org/oklein/dune-randomfield.git
#    git clone https://parcomp-git.iwr.uni-heidelberg.de/Teaching/dune-funcep.git
#    git clone https://parcomp-git.iwr.uni-heidelberg.de/Teaching/dune-npde.git
#    git clone https://parcomp-git.iwr.uni-heidelberg.de/Teaching/dune-parsolve.git
#    git clone https://parcomp-git.iwr.uni-heidelberg.de/peter/dune-padua.git
#    git clone https://parcomp-git.iwr.uni-heidelberg.de/peter/dune-dr.git
    git clone https://parcomp-git.iwr.uni-heidelberg.de/peter/dune-hydro.git
#    git clone https://parcomp-git.iwr.uni-heidelberg.de/peter/dune-ftworth.git
#    git clone https://gitlab.dune-project.org/dominic/dune-fiber-elasticity.git
    git clone https://gitlab.dune-project.org/oklein/dune-nonlinopt.git
    
    # install whitespace hook
    ./dune-common/bin/dunecontrol vcsetup
}

# generate options file for Dune build system
# generate options file for Dune build system
generate_optsfile ()
{
    CC=$(which $CC)
    CXX=$(which $CXX)
    F77=$(which $F77)
echo "CMAKE_FLAGS=\"
-DCMAKE_C_COMPILER='$CC'
-DCMAKE_CXX_COMPILER='$CXX'
-DCMAKE_Fortran_COMPILER='$F77'
-DCMAKE_CXX_FLAGS_RELEASE='-O3 -DNDEBUG -g0 -funroll-loops -ftemplate-depth=5120 -march=native -Wa,-q'
-DCMAKE_BUILD_TYPE=Release
-DDUNE_SYMLINK_TO_SOURCE_TREE=1
-DARPACK_ROOT=$EXTERNAL_HOME/arpack-ng
-DARPACKPP_INCLUDE_DIR=$EXTERNAL_HOME/arpackpp/include/arpackpp/
\"" > release.opts
}

extra_flags ()
{
# These are just more flags that you may put in generate_optsfile ()
echo "CMAKE_FLAGS=\"
-DVCL_ROOT=$EXTERNAL_HOME/vcl
-DHDF5_ROOT=$EXTERNAL_HOME/hdf5
-DFFTW3_ROOT_DIR=$EXTERNAL_HOME/fftw3
-DCMAKE_DISABLE_FIND_PACKAGE_ParMETIS=1
-DCMAKE_DISABLE_FIND_PACKAGE_SuperLU=1
\""
}

# create build script
generate_buildscript ()
{
echo '#!/bin/bash
ROOT=$(pwd)
BUILDDIR=$ROOT/release-build
OPTSFILE=release.opts
MODULES="common geometry uggrid grid alugrid istl localfunctions typetree functions testtools pdelab pdelab-tutorials"
for i in $MODULES; do
    echo build $i
    ./dune-common/bin/dunecontrol --builddir=$BUILDDIR  --opts=$OPTSFILE --only=dune-$i all
done
' > buildmodules.sh
chmod 755 buildmodules.sh
}

# now run your selection of commands
#install_vcl
install_arpack-ng
#install_fftw3
#install_hdf5
checkout_dune_core
checkout_dune_pdelab
checkout_dune_pdelab_applications
generate_optsfile
generate_buildscript
