# Pull a debian image
#
#    docker pull debian:stable
#
# Start it
#
#    docker run -it --name debian debian:stable
#
# Inside the container, run the following commands to install firedrkae.
#
# Subsequently, you can start and attach to this container
#
#    docker start debian
#    docker attach debian
#
# Activate the env
#
#    . /root/firedrake/bin/activate
#
# Add to your .bashrc file
# export OMP_NUM_THREADS=1
#
# CREATE AN IMAGE
#
# You can create an image out of this container
#
#    docker commit debian firedrake
#
# Run the image (you may want to share some host folder)
#
#    docker run -it --name firedrake -p 8888:8888 firedrake:latest
#
# You can start jupyter inside container 
#
#    . /root/firedrake/bin/activate
#    jupyter-lab --ip 0.0.0.0 --port 8888 --no-browser --allow-root
#
# and access it from the host at http://localhost:8888/tree
cd
apt update
apt install curl git python3 python3-venv vim
curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/refs/tags/2025.4.2/scripts/firedrake-configure
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
