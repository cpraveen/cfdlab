echo export PETSC_DIR=$(spack location -i petsc)
echo export SLEPC_DIR=$(spack location -i slepc)
echo export HDF5_DIR=$(spack location -i hdf5)
echo export METIS_DIR=$(spack location -i metis)

PACKAGES="cmake openmpi "

for pkg in $PACKAGES
do
   echo PATH=\$PATH:$(spack location -i $pkg)/bin
done
