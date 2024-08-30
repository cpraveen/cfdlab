unset PETSC_DIR
unset SLEPC_DIR
unset HDF5_DIR

PATH=$PATH:`spack location -i autoconf`/bin
PATH=$PATH:`spack location -i automake`/bin
PATH=$PATH:`spack location -i gmake`/bin
PATH=$PATH:`spack location -i libtool`/bin
PATH=$PATH:`spack location -i ninja`/bin
echo $PATH

#PETSC_CONFIGURE_OPTIONS="--with-hdf5-dir=$HDF5_DIR"
#
export PETSC_CONFIGURE_OPTIONS="--download-zlib"

python ./firedrake-install --show-dependencies

python ./firedrake-install --venv-name firedrake \
                           --slepc \
                           --with-parmetis \
                           --mpicc mpicc \
                           --mpicxx mpicxx \
                           --mpif90 mpif90 \
                           --mpiexec mpiexec \
                           --with-blas /opt/spack \
                           --no-package-manager
