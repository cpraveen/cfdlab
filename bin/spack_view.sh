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

#read -p "Link gcc (y/n) ? " CONT
#if [ "$CONT" = "y" ]; then
#   dolink gcc
#   echo "Deleting c++, cpp from spack view"
#   rm -f $SPACK_VIEW/bin/c++
#   rm -f $SPACK_VIEW/bin/cpp
#
#   read -p "Continue (y/n) ? " CONT
#   if [ "$CONT" = "n" ]; then
#      exit
#   fi
#fi

dolink adol-c
dolink arborx
dolink arpack-ng
dolink assimp
dolink boost
dolink cmake
dolink eigen
dolink ginkgo
dolink gmake
dolink gsl
dolink hdf5
dolink kokkos
dolink metis
dolink mpich
dolink muparser
dolink mumps
dolink netlib-scalapack
dolink oce
dolink openblas
dolink openmpi
dolink parmetis
dolink p4est
dolink petsc
dolink slepc
dolink suite-sparse
dolink sundials
dolink symengine
dolink taskflow
dolink tbb
dolink trilinos
dolink zlib
dolink zlib-ng
