# This script runs cmake for deal.II in combination with Spack
#
# Extract deal.II, change into the deal.II directory and
#     mkdir build && cd build
# Then run this script inside "build" directory
#     sh <path to this script> <path to spack view> <path to deal.II install dir>
# Then compile deal.II (set -j# based on cores you have for parallel compile)
#     make -j8 all
# Install deal.II
#     make install
#

PN=`basename "$0"`

usage () {
   echo >&2 "usage: $PN <spack view dir> <deal.II install dir>

   Configure deal.II using dependencies installed via Spack view
   example: sh $PN /opt/spack /home/praveen/Applications/deal.II/install"

    exit 1
}

if [ $# -lt 2 ]
then
	usage
fi

# Location of Spack
SPACK_VIEW_DIR=$1

# Where do you want to install deal.II
DEAL_II_DIR=$2


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

MPI_DIR=$SPACK_VIEW_DIR

# Get c/c++/fortran from mpi wrappers
CC=`$MPI_DIR/bin/mpicc -show   | awk '{print $1}'`
CXX=`$MPI_DIR/bin/mpic++ -show | awk '{print $1}'`
FC=`$MPI_DIR/bin/mpif90 -show  | awk '{print $1}'`

BLAS_DIR=$SPACK_VIEW_DIR

# Create makefile
$SPACK_VIEW_DIR/bin/cmake  \
-DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR \
-DCMAKE_FIND_FRAMEWORK=LAST  \
-DCMAKE_BUILD_TYPE=DebugRelease  \
-DDEAL_II_COMPONENT_EXAMPLES=ON  \
-DDEAL_II_COMPILE_EXAMPLES=OFF \
-DDEAL_II_COMPONENT_DOCUMENTATION=OFF  \
-DDEAL_II_WITH_LAPACK=ON \
-DLAPACK_DIR=$BLAS_DIR  \
-DBLAS_DIR=$BLAS_DIR  \
-DDEAL_II_WITH_BOOST=ON \
-DBOOST_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_ARBORX=ON \
-DARBORX_DIR=$SPACK_VIEW_DIR \
-DMUPARSER_DIR=$SPACK_VIEW_DIR \
-DUMFPACK_DIR=$SPACK_VIEW_DIR \
-DTBB_DIR=$SPACK_VIEW_DIR \
-DZLIB_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_MPI=ON  \
-DCMAKE_C_COMPILER=$CC \
-DCMAKE_CXX_COMPILER=$CXX \
-DCMAKE_Fortran_COMPILER=$FC \
-DMPI_C_COMPILER=$MPI_DIR/bin/mpicc  \
-DMPI_CXX_COMPILER=$MPI_DIR/bin/mpic++  \
-DMPI_Fortran_COMPILER=$MPI_DIR/bin/mpif90  \
-DDEAL_II_CXX_FLAGS="-march=native -mtune=native" \
-DDEAL_II_CXX_FLAGS_RELEASE="-O3" \
-DGSL_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_GSL=ON  \
-DGMP_INCLUDE_DIR=$SPACK_VIEW_DIR \
-DGMP_LIBRARIES=$SPACK_VIEW_DIR \
-DMPFR_INCLUDE_DIR=$SPACK_VIEW_DIR \
-DMPFR_LIBRARIES=$SPACK_VIEW_DIR \
-DCGAL_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_CGAL=ON  \
-DHDF5_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_HDF5=ON  \
-DP4EST_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_P4EST=ON  \
-DPETSC_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_PETSC=ON  \
-DSLEPC_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_SLEPC=ON  \
-DTRILINOS_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_TRILINOS=ON  \
-DMETIS_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_METIS=ON  \
-DMUMPS_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_MUMPS=ON  \
-DARPACK_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_ARPACK=ON  \
-DDEAL_II_ARPACK_WITH_PARPACK=ON  \
-DOPENCASCADE_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_OPENCASCADE=ON  \
-DGMSH_DIR=$GMSH_DIR \
-DDEAL_II_WITH_GMSH=$DEAL_II_WITH_GMSH  \
-DASSIMP_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_ASSIMP=ON  \
-DSUNDIALS_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_SUNDIALS=ON  \
-DADOLC_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_ADOLC=ON  \
-DSCALAPACK_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_SCALAPACK=ON  \
-DDEAL_II_WITH_GINKGO=ON  \
-DGINKGO_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_SYMENGINE=ON  \
-DSYMENGINE_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_TASKFLOW=ON  \
-DTASKFLOW_DIR=$SPACK_VIEW_DIR \
-DDEAL_II_WITH_VTK=OFF  \
../

echo "Check detailed.log and CMakeCache.txt"
echo "*** Add this to your profile ***"
echo "   export DEAL_II_DIR=$DEAL_II_DIR"
