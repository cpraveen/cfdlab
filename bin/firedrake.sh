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
python3 ./firedrake-install --show-dependencies
#python3 ./firedrake-install --show-dependencies --honour-petsc-dir

python3 ./firedrake-install --venv-name firedrake \
                            --slepc \
                            --with-parmetis \
                            --mpicc /opt/spack/bin/mpicc \
                            --mpicxx /opt/spack/bin/mpicxx \
                            --mpif90 /opt/spack/bin/mpif90 \
                            --mpiexec /opt/spack/bin/mpiexec \
                            --mpihome /opt/spack \
                            --with-blas /opt/spack \
                            --no-package-manager
