# See
#    https://cpraveen.github.io/comp/firedrake
#
# Exit on error
set -e
rm -rf /root/.cache
apt update
apt install curl git python3 python3-venv vim
cd /root
curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/refs/tags/2025.4.3/scripts/firedrake-configure
apt install $(python3 firedrake-configure --show-system-packages)
git clone --branch $(python3 firedrake-configure --show-petsc-version) https://gitlab.com/petsc/petsc.git
cd petsc
python3 ../firedrake-configure --show-petsc-configure-options | xargs -L1 ./configure
make PETSC_DIR=/root/petsc PETSC_ARCH=arch-firedrake-default all
make PETSC_DIR=/root/petsc PETSC_ARCH=arch-firedrake-default check
cd ..
python3 -m venv firedrake
. firedrake/bin/activate
pip cache purge
pip cache remove mpi4py
pip cache remove petsc4py
pip cache remove h5py
pip cache remove slepc4py
pip cache remove libsupermesh
pip cache remove firedrake
export $(python3 firedrake-configure --show-env)
pip install --no-binary h5py 'firedrake[check]'
pip install numpy scipy sympy matplotlib jupyterlab vtk
firedrake-check

# Do some cleanup

firedrake-clean
apt autoremove --purge
apt clean
rm -rf /var/lib/apt/lists/*
rm -rf `find /root/petsc -name .git`
rm -rf /root/petsc/src/docs
rm -f  /root/petsc/src/**/tutorials/output/*
rm -f  /root/petsc/src/**/tests/output/*
