#!/bin/bash
ROOT=$(pwd)
BUILDDIR=$ROOT/release-build
OPTSFILE=release.opts
MODULES="common geometry uggrid grid alugrid spgrid multidomaingrid istl localfunctions typetree functions testtools pdelab pdelab-tutorials fem ftworth moes"
for i in $MODULES; do
    echo build $i
    ./dune-common/bin/dunecontrol --builddir=$BUILDDIR  --opts=$OPTSFILE --only=dune-$i all
done

