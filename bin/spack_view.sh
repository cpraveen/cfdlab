#!/bin/bash

PN=`basename "$0"`

usage () {
   echo >&2 "usage: $PN <spack view dir>

   Install spack view into specified directory.
   example: sh $PN /opt/spack"

    exit 1
}

if [ $# -lt 1 ]
then
	usage
fi

# Set this to your actual spack view dir
SPACK_VIEW=$1

COMMAND="spack view -v -d no"
ACTION=symlink

dolink()
{
   $COMMAND $ACTION $SPACK_VIEW $1
   # Delete some useless stuff
   rm -f $SPACK_VIEW/include/index.html
   rm -f $SPACK_VIEW/share/info/dir
}

if [ ! -d "$SPACK_VIEW" ]; then
   echo "Directory $SPACK_VIEW does not exist"
   exit
fi

# Delete all existing stuff
rm -rf $SPACK_VIEW/*
rm -rf $SPACK_VIEW/.spack
rm -rf $SPACK_VIEW/.spack-empty
rm -rf $SPACK_VIEW/.nagged

dolink adol-c
dolink arpack-ng
dolink assimp
dolink boost
dolink bzip2
dolink cmake
dolink ginkgo
dolink gmsh
dolink gsl
dolink hdf5
dolink metis
dolink mpich
dolink muparser
dolink mumps
dolink nanoflann
dolink netcdf-c
dolink netcdf-cxx
dolink netlib-scalapack
dolink oce
dolink openblas
dolink openmpi
dolink p4est
dolink petsc
dolink slepc
dolink suite-sparse
dolink sundials
dolink symengine
dolink tbb
dolink trilinos
dolink zlib
