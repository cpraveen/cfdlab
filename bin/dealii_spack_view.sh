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

# Determine os
OS=`uname -s`
echo "OS is $OS"
if [ $OS = "Darwin" ]; then
   LAPACK_LIBRARY=libopenblas.dylib
elif [ $OS = "Linux" ]; then
   LAPACK_LIBRARY=libopenblas.so
else
   echo "Unable to determine LAPACK_LIBRARY"
   exit
fi

${SPACK_VIEW_DIR}/bin/cmake  \
-DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR \
-DCMAKE_FIND_FRAMEWORK=LAST  \
-DCMAKE_BUILD_TYPE=DebugRelease  \
-DDEAL_II_COMPONENT_EXAMPLES=ON  \
-DDEAL_II_COMPILE_EXAMPLES=OFF \
-DDEAL_II_WITH_ARBORX=ON  \
-DDEAL_II_WITH_LAPACK=ON \
-DLAPACK_INCLUDE_DIRS=${SPACK_VIEW_DIR}/include  \
-DLAPACK_LIBRARIES=${SPACK_VIEW_DIR}/lib/${LAPACK_LIBRARY}  \
-DARBORX_DIR=${SPACK_VIEW_DIR}  \
-DBOOST_DIR=${SPACK_VIEW_DIR}  \
-DMUPARSER_DIR=${SPACK_VIEW_DIR}  \
-DUMFPACK_DIR=${SPACK_VIEW_DIR}  \
-DTBB_DIR=${SPACK_VIEW_DIR}  \
-DZLIB_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_MPI:BOOL=ON  \
-DCMAKE_C_COMPILER=${SPACK_VIEW_DIR}/bin/mpicc  \
-DCMAKE_CXX_COMPILER=${SPACK_VIEW_DIR}/bin/mpic++  \
-DCMAKE_Fortran_COMPILER=${SPACK_VIEW_DIR}/bin/mpif90  \
-DDEAL_II_CXX_FLAGS="-march=native -std=c++17" \
-DDEAL_II_CXX_FLAGS_RELEASE="-O3" \
-DGSL_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_GSL:BOOL=ON  \
-DHDF5_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_HDF5:BOOL=ON  \
-DP4EST_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_P4EST:BOOL=ON  \
-DPETSC_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_PETSC:BOOL=ON  \
-DSLEPC_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_SLEPC:BOOL=ON  \
-DTRILINOS_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_TRILINOS:BOOL=ON  \
-DMETIS_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_METIS:BOOL=ON  \
-DDEAL_II_COMPONENT_DOCUMENTATION=OFF  \
-DARPACK_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_ARPACK=ON  \
-DDEAL_II_ARPACK_WITH_PARPACK=ON  \
-DOPENCASCADE_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_OPENCASCADE=ON  \
-DGMSH_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_GMSH=ON  \
-DASSIMP_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_ASSIMP=ON  \
-DSUNDIALS_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_SUNDIALS=ON  \
-DADOLC_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_ADOLC=ON  \
-DSCALAPACK_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_SCALAPACK=ON  \
-DDEAL_II_WITH_GINKGO=ON  \
-DGINKGO_DIR=${SPACK_VIEW_DIR}  \
-DDEAL_II_WITH_SYMENGINE=ON  \
-DSYMENGINE_DIR=${SPACK_VIEW_DIR}  \
../

echo "*** Add this to your profile ***"
echo "export DEAL_II_DIR=$DEAL_II_DIR"
