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

   Configure deal.II using dependencies installed via Spack view
   example: $PN /home/praveen/Applications/deal.II/install"

    exit 1
}

if [ $# -lt 1 ]
then
	usage
fi

# Where do you want to install deal.II
DEAL_II_DIR=$1

# Determine os
OS=`uname -s`
echo "OS is $OS"
if [ $OS = "Darwin" ]; then
   LAPACK_LIBRARY=libopenblas.dylib
elif [ $OS = "Linux" ]; then
   LAPACK_LIBRARY=libopenblas.so
else
   echo "Problem with LAPACK_LIBRARY"
   exit
fi


$(spack location -i cmake)/bin/cmake  \
-DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR \
-DCMAKE_FIND_FRAMEWORK=LAST  \
-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON  \
-DCMAKE_BUILD_TYPE=DebugRelease  \
-DDEAL_II_COMPONENT_EXAMPLES=ON  \
-DDEAL_II_COMPILE_EXAMPLES=OFF \
-DDEAL_II_WITH_ARBORX=ON \
-DDEAL_II_WITH_LAPACK=ON \
-DLAPACK_INCLUDE_DIRS=$(spack location -i openblas)/include  \
-DLAPACK_LIBRARIES=$(spack location -i openblas)/lib/${LAPACK_LIBRARY}  \
-DBOOST_DIR=$(spack location -i boost)  \
-DARBORX_DIR=$(spack location -i arborx)  \
-DMUPARSER_DIR=$(spack location -i muparser)  \
-DUMFPACK_DIR=$(spack location -i suite-sparse)  \
-DTBB_DIR=$(spack location -i intel-tbb)  \
-DZLIB_DIR=$(spack location -i zlib)  \
-DDEAL_II_WITH_MPI:BOOL=ON  \
-DCMAKE_C_COMPILER=$(spack location -i openmpi)/bin/mpicc  \
-DCMAKE_CXX_COMPILER=$(spack location -i openmpi)/bin/mpic++  \
-DCMAKE_Fortran_COMPILER=$(spack location -i openmpi)/bin/mpif90  \
-DDEAL_II_CXX_FLAGS="-march=native -std=c++17" \
-DDEAL_II_CXX_FLAGS_RELEASE="-O3" \
-DGSL_DIR=$(spack location -i gsl)  \
-DDEAL_II_WITH_GSL:BOOL=ON  \
-DHDF5_DIR=$(spack location -i hdf5)  \
-DDEAL_II_WITH_HDF5:BOOL=ON  \
-DP4EST_DIR=$(spack location -i p4est)  \
-DDEAL_II_WITH_P4EST:BOOL=ON  \
-DPETSC_DIR=$(spack location -i petsc)  \
-DDEAL_II_WITH_PETSC:BOOL=ON  \
-DSLEPC_DIR=$(spack location -i slepc)  \
-DDEAL_II_WITH_SLEPC:BOOL=ON  \
-DTRILINOS_DIR=$(spack location -i trilinos)  \
-DDEAL_II_WITH_TRILINOS:BOOL=ON  \
-DMETIS_DIR=$(spack location -i metis)  \
-DDEAL_II_WITH_METIS:BOOL=ON  \
-DDEAL_II_COMPONENT_DOCUMENTATION=OFF  \
-DARPACK_DIR=$(spack location -i arpack-ng)  \
-DDEAL_II_WITH_ARPACK=ON  \
-DDEAL_II_ARPACK_WITH_PARPACK=ON  \
-DOPENCASCADE_DIR=$(spack location -i oce)  \
-DDEAL_II_WITH_OPENCASCADE=ON  \
-DGMSH_DIR=$(spack location -i gmsh)  \
-DDEAL_II_WITH_GMSH=ON  \
-DASSIMP_DIR=$(spack location -i assimp)  \
-DDEAL_II_WITH_ASSIMP=ON  \
-DSUNDIALS_DIR=$(spack location -i sundials)  \
-DDEAL_II_WITH_SUNDIALS=ON  \
-DADOLC_DIR=$(spack location -i adol-c)  \
-DDEAL_II_WITH_ADOLC=ON  \
-DSCALAPACK_DIR=$(spack location -i netlib-scalapack)  \
-DDEAL_II_WITH_SCALAPACK=ON  \
-DDEAL_II_WITH_GINKGO=ON  \
-DGINKGO_DIR=$(spack location -i ginkgo)  \
-DDEAL_II_WITH_SYMENGINE=ON  \
-DSYMENGINE_DIR=$(spack location -i symengine)  \
../

echo "*** Add this to your profile ***"
echo "export DEAL_II_DIR=$DEAL_II_DIR"
