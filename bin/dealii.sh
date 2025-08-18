# This script runs cmake for deal.II in combination with Spack
#
# Extract deal.II, change into the deal.II directory and
#     mkdir build && cd build
# Then run this script inside "build" directory
#     sh <path to this script> <path to deal.II install dir>
# Then compile deal.II (set -j# based on cores you have for parallel compile)
#     make -j8 all
# Install deal.II
#     make install
#

PN=`basename "$0"`

usage () {
   echo >&2 "usage: $PN <deal.II install dir>

   Configure deal.II assuming everything is in path.
   example: $PN /home/praveen/Applications/deal.II/install"

    exit 1
}

if [ $# -ne 1 ]
then
	usage
fi

# Where do you want to install deal.II
DEAL_II_DIR=$1

# Determine os
OS=`uname -s`
echo "OS is $OS"
if [ $OS = "Darwin" ]; then
   export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
fi

# We dont want to install gmsh with spack since it has too many dependencies.
# We want to use externally installed gmsh, so check here.
DEAL_II_WITH_GMSH=ON
if [ -z "$GMSH_DIR" ]; then
   echo "GMSH_DIR is not set"
   read -p "Do you want to continue (y/n) ? " CONT
   if [ "$CONT" = "n" ]; then
      exit
   fi
   DEAL_II_WITH_GMSH=OFF
   GMSH_DIR=/tmp
fi

# Get c/c++/fortran from mpi wrappers
CC=`mpicc -show   | awk '{print $1}'`
CXX=`mpic++ -show | awk '{print $1}'`
FC=`mpif90 -show  | awk '{print $1}'`

# Create makefile
cmake  \
-DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR \
-DCMAKE_FIND_FRAMEWORK=LAST  \
-DCMAKE_BUILD_TYPE=DebugRelease  \
-DDEAL_II_COMPONENT_EXAMPLES=ON  \
-DDEAL_II_COMPILE_EXAMPLES=OFF \
-DDEAL_II_COMPONENT_DOCUMENTATION=OFF  \
-DDEAL_II_WITH_LAPACK=ON \
-DDEAL_II_WITH_BOOST=ON \
-DDEAL_II_WITH_ARBORX=ON \
-DDEAL_II_WITH_MPI=ON  \
-DCMAKE_C_COMPILER=$CC \
-DCMAKE_CXX_COMPILER=$CXX \
-DCMAKE_Fortran_COMPILER=$FC \
-DMPI_C_COMPILER=mpicc  \
-DMPI_CXX_COMPILER=mpic++  \
-DMPI_Fortran_COMPILER=mpif90  \
-DDEAL_II_CXX_FLAGS="-march=native -mtune=native" \
-DDEAL_II_CXX_FLAGS_RELEASE="-O3" \
-DDEAL_II_WITH_GSL=ON  \
-DDEAL_II_WITH_CGAL=ON  \
-DDEAL_II_WITH_HDF5=ON  \
-DDEAL_II_WITH_P4EST=ON  \
-DDEAL_II_WITH_PETSC=ON  \
-DDEAL_II_WITH_SLEPC=ON  \
-DDEAL_II_WITH_TRILINOS=ON  \
-DDEAL_II_WITH_METIS=ON  \
-DDEAL_II_WITH_MUMPS=ON  \
-DDEAL_II_WITH_ARPACK=ON  \
-DDEAL_II_ARPACK_WITH_PARPACK=ON  \
-DDEAL_II_WITH_OPENCASCADE=ON  \
-DGMSH_DIR=$GMSH_DIR \
-DDEAL_II_WITH_GMSH=$DEAL_II_WITH_GMSH  \
-DDEAL_II_WITH_ASSIMP=ON  \
-DDEAL_II_WITH_SUNDIALS=ON  \
-DDEAL_II_WITH_ADOLC=ON  \
-DDEAL_II_WITH_SCALAPACK=ON  \
-DDEAL_II_WITH_GINKGO=ON  \
-DDEAL_II_WITH_SYMENGINE=ON  \
-DDEAL_II_WITH_TASKFLOW=ON  \
-DDEAL_II_WITH_VTK=OFF  \
../

echo "Check detailed.log and CMakeCache.txt for paths you dont want to include."
echo "E.g., if you have miniforge/anaconda, check they are not used."
echo "Make sure mpicc, mpirun, etc. are present in PATH."
echo "*** Add this to your profile ***"
echo "   export DEAL_II_DIR=$DEAL_II_DIR"
