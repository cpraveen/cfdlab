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
   LAPACK_LIBRARY=libopenblas.dylib
elif [ $OS = "Linux" ]; then
   LAPACK_LIBRARY=libopenblas.so
else
   echo "Problem with LAPACK_LIBRARY"
   exit
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

# Select mpich or openmpi
OPENMPI=`spack location -i openmpi 2> /dev/null`
MPICH=`spack location -i mpich 2> /dev/null`

if [ $OPENMPI ]; then
   MPI_DIR=$OPENMPI
else
   MPI_DIR=$MPICH
fi

#read -p "mpich or openmpi ? " MPI
#if [ $MPI = "mpich" ]; then
#  MPI_DIR=`spack location -i mpich`
#elif [ $MPI = "openmpi" ]; then
#  MPI_DIR=`spack location -i openmpi`
#else
#   echo "Unknown MPI specified"
#   exit
#fi

if [ -d "$MPI_DIR" ]; then
   echo "MPI_DIR = " $MPI_DIR
else
   echo "MPI_DIR not found"
   exit
fi

`spack location -i cmake`/bin/cmake  \
-DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR \
-DCMAKE_FIND_FRAMEWORK=LAST  \
-DCMAKE_BUILD_TYPE=DebugRelease  \
-DDEAL_II_COMPONENT_EXAMPLES=ON  \
-DDEAL_II_COMPILE_EXAMPLES=OFF \
-DDEAL_II_WITH_LAPACK=ON \
-DLAPACK_INCLUDE_DIRS=`spack location -i openblas`/include  \
-DLAPACK_LIBRARIES=`spack location -i openblas`/lib/${LAPACK_LIBRARY}  \
-DBOOST_DIR=`spack location -i boost`  \
-DDEAL_II_WITH_ARBORX=ON \
-DARBORX_DIR=`spack location -i arborx`  \
-DMUPARSER_DIR=`spack location -i muparser`  \
-DUMFPACK_DIR=`spack location -i suite-sparse`  \
-DTBB_DIR=`spack location -i intel-tbb`  \
-DZLIB_DIR=`spack location -i zlib-ng`  \
-DDEAL_II_WITH_MPI=ON  \
-DCMAKE_C_COMPILER=$MPI_DIR/bin/mpicc  \
-DCMAKE_CXX_COMPILER=$MPI_DIR/bin/mpic++  \
-DCMAKE_Fortran_COMPILER=$MPI_DIR/bin/mpif90  \
-DDEAL_II_CXX_FLAGS="-march=native -mtune=native -std=c++17" \
-DDEAL_II_CXX_FLAGS_RELEASE="-O3" \
-DGSL_DIR=`spack location -i gsl`  \
-DDEAL_II_WITH_GSL=ON  \
-DGMP_INCLUDE_DIR=`spack location -i gmp` \
-DGMP_LIBRARIES=`spack location -i gmp` \
-DMPFR_INCLUDE_DIR=`spack location -i mpfr` \
-DMPFR_LIBRARIES=`spack location -i mpfr` \
-DCGAL_DIR=`spack location -i cgal`  \
-DDEAL_II_WITH_CGAL=ON  \
-DHDF5_DIR=`spack location -i hdf5`  \
-DDEAL_II_WITH_HDF5=ON  \
-DP4EST_DIR=`spack location -i p4est`  \
-DDEAL_II_WITH_P4EST=ON  \
-DPETSC_DIR=`spack location -i petsc`  \
-DDEAL_II_WITH_PETSC=ON  \
-DSLEPC_DIR=`spack location -i slepc`  \
-DDEAL_II_WITH_SLEPC=ON  \
-DTRILINOS_DIR=`spack location -i trilinos`  \
-DDEAL_II_WITH_TRILINOS=ON  \
-DMETIS_DIR=`spack location -i metis`  \
-DDEAL_II_WITH_METIS=ON  \
-DDEAL_II_COMPONENT_DOCUMENTATION=OFF  \
-DARPACK_DIR=`spack location -i arpack-ng`  \
-DDEAL_II_WITH_ARPACK=ON  \
-DDEAL_II_ARPACK_WITH_PARPACK=ON  \
-DOPENCASCADE_DIR=`spack location -i oce`  \
-DDEAL_II_WITH_OPENCASCADE=ON  \
-DGMSH_DIR=$GMSH_DIR \
-DDEAL_II_WITH_GMSH=$DEAL_II_WITH_GMSH  \
-DASSIMP_DIR=`spack location -i assimp`  \
-DDEAL_II_WITH_ASSIMP=ON  \
-DSUNDIALS_DIR=`spack location -i sundials`  \
-DDEAL_II_WITH_SUNDIALS=ON  \
-DADOLC_DIR=`spack location -i adol-c`  \
-DDEAL_II_WITH_ADOLC=ON  \
-DSCALAPACK_DIR=`spack location -i netlib-scalapack`  \
-DDEAL_II_WITH_SCALAPACK=ON  \
-DDEAL_II_WITH_GINKGO=ON  \
-DGINKGO_DIR=`spack location -i ginkgo`  \
-DDEAL_II_WITH_SYMENGINE=ON  \
-DSYMENGINE_DIR=`spack location -i symengine`  \
-DDEAL_II_WITH_TASKFLOW=ON  \
-DTASKFLOW_DIR=`spack location -i taskflow`  \
-DDEAL_II_WITH_VTK=OFF  \
../

echo "*** Add this to your profile ***"
echo "export DEAL_II_DIR=$DEAL_II_DIR"
