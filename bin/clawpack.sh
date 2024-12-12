# Install clawpack from source, uses conda
# Tested to work with bash and zsh

# Exit on error
set -e

export GIT_ADVICE=0

# Check for git
if [ ! type  git &> /dev/null ]; then
   echo "git is not found, install/add it to your path and try again"
   exit
fi

# Check for conda
if [ ! type  conda &> /dev/null ]; then
   echo "conda is not found, install/add it to your path and try again"
   exit
fi

# Check that CLAW is set
if [ -z $CLAW ]; then
   echo "Set CLAW to full path of clawpack directory"
   echo "E.g., export CLAW=$HOME/Applications/clawpack"
   echo "Git sources of clawpack will be cloned into this."
   echo "If this directory exists, it will be deleted."
   exit
fi

# Check $CLAW is writable
if [ ! -w "$(dirname "$CLAW")" ]; then
   echo "CLAW = $CLAW" 
   echo "No permission to create/write into CLAW."
   echo "Set CLAW to some location with write permission."
   exit
fi

# Need to specify clawpack version to install
if [ $# -eq 0 ]; then
    echo "Clawpack version is not supplied"
    echo "Example: run this script like this"
    echo "   bash ./clawpack.sh v5.9.2"
    exit
fi

# clawpack version
VERSION=$1

# Name of conda environment, you can change this if you want
ENV=claw

echo "Installing in conda env: $ENV"
echo "CLAW = $CLAW"
echo "Clawpack sources will be downloaded to above directory."
if [ -d $CLAW ]; then
   GIT=`cd $CLAW && git describe --all 2> /dev/null`
   echo "WARNING: Directory"
   echo "            $CLAW"
   echo "         exists\c"
   if [ $GIT ]; then
      echo " with Clawpack git version: $GIT"
   else
      echo " but there is no Clawpack git version."
   fi
   read -p "Do you want to delete and checkout again ? (y/n/ctr-c) " update_claw
else
   echo "Directory $CLAW will be created"
   update_claw="y"
fi

if [[ "$update_claw" == "y" ]]; then

echo "----------------------------------------------------------------------"
echo "Checking out clawpack source from git"
echo "----------------------------------------------------------------------"
rm -rf $CLAW
git clone https://github.com/clawpack/clawpack.git $CLAW
cd $CLAW
git clone --recursive https://github.com/clawpack/apps
git checkout $VERSION
git submodule init
git submodule update

else

echo "----------------------------------------------------------------------"
echo "Checking out clawpack version " $VERSION
echo "----------------------------------------------------------------------"
cd $CLAW
git checkout master
git pull
git checkout $VERSION
git submodule init
git submodule update

fi

# Install needed packages inside conda env
eval "$(conda shell.bash hook)"

find_in_conda_env(){
    conda env list | grep "${@} " >/dev/null 2>/dev/null
}

if find_in_conda_env $ENV ; then
   echo "----------------------------------------------------------------------"
   echo "Conda env $ENV exists"
   echo "----------------------------------------------------------------------"
   read -p "Delete existing conda env ? (y/n/ctr-c) " use_env
   if [[ "$use_env" == "n" ]]; then
      conda activate $ENV
   else
      conda remove -y -n claw --all
      conda create -y -n $ENV
      conda activate $ENV
   fi
else
   echo "----------------------------------------------------------------------"
   echo "Creating conda env $ENV and installing packages"
   echo "----------------------------------------------------------------------"
   conda create -y -n $ENV
   conda activate $ENV
fi

PACKAGES="gfortran ipython ipywidgets jupyterlab matplotlib meson-python \
          ninja nose numpy petsc4py pip pkg-config pytest \
          scipy seaborn six spin sympy"

# If gfortran is not found, then install it
#if command -v  gfortran &> /dev/null
#then
#   echo "gfortran already exists: " `command -v gfortran`
#   gfortran --version | head -n1
#   read -p "Install gfortran with conda ? (y/n/ctr-c) " fortran
#   if [[ "$fortran" == "y" ]]; then
#      PACKAGES="$PACKAGES gfortran"
#   fi
#else
#   echo "gfortran not found, will be installed in the conda env"
#   PACKAGES="$PACKAGES gfortran"
#fi

conda install -y -c conda-forge $PACKAGES

export FC=$CONDA_PREFIX/bin/gfortran
if [[ $OSTYPE == 'darwin'* ]]; then
   export CC=$CONDA_PREFIX/bin/clang
else
   export CC=$CONDA_PREFIX/bin/cc
fi

# Build clawpack
echo "----------------------------------------------------------------------"
echo "Building clawpack"
echo "----------------------------------------------------------------------"
cd $CLAW
rm -rf build
pip install --user --no-build-isolation -e .

echo "----------------------------------------------------------------------"
echo "Installation successful, try to run a test"
echo "   conda activate claw"
echo "   cd /tmp"
echo "   python $CLAW/pyclaw/examples/advection_1d/advection_1d.py iplot=1"
echo "It is better not to run inside $CLAW or to modify the code there."
echo "----------------------------------------------------------------------"
