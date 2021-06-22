echo export PETSC_DIR=$(spack location -i petsc)
echo export SLEPC_DIR=$(spack location -i slepc)
echo export HDF5_DIR=$(spack location -i hdf5)
echo export METIS_DIR=$(spack location -i metis)
echo export EIGEN_DIR=$(spack location -i eigen)
echo export BLAS_DIR=$(spack location -i openblas)
echo export LAPACK_DIR=$(spack location -i openblas)

PACKAGES="cmake openmpi "

for pkg in $PACKAGES
do
   echo PATH=\$PATH:$(spack location -i $pkg)/bin
done
