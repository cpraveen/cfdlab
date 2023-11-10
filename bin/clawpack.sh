# Install clawpack from source, uses conda
set -e

if [ -z `which git` ]; then
   echo "git is not found, install/add it to your path and try again"
   exit
fi

if [ -z `which conda` ]; then
   echo "conda is not found, install/add it to your path and try again"
   exit
fi

if [ -z $CLAW ]; then
   echo "Set CLAW to full path of clawpack directory"
   echo "E.g., export CLAW=$HOME/Applications/clawpack"
   echo "Git sources of clawpack will be cloned into this."
   echo "If this directory exists, it will be deleted."
   exit
fi

if [ $# -eq 0 ]
  then
    echo "Clawpack version is not supplied"
    echo "Example: run this script like this"
    echo "   sh ./clawpack.sh v5.9.2"
    exit
fi

VERSION=$1

# Name of conda environment
ENV=claw

echo "Will install using conda env: $ENV"
echo "Directory $CLAW will be deleted if it exists."
read -p "Press enter to continue or control-c to quit "

echo "----------------------------------------------------------------------"
echo "Checking out clawpack source"
echo "----------------------------------------------------------------------"
rm -rf $CLAW
git clone https://github.com/clawpack/clawpack.git $CLAW
cd $CLAW
git clone --recursive https://github.com/clawpack/apps
git checkout $VERSION
git submodule init
git submodule update

eval "$(conda shell.bash hook)"

find_in_conda_env(){
    conda env list | grep "${@}" >/dev/null 2>/dev/null
}

if find_in_conda_env $ENV ; then
   echo "----------------------------------------------------------------------"
   echo "Conda env $ENV exists, activating"
   echo "----------------------------------------------------------------------"
   conda activate $ENV
else
   echo "----------------------------------------------------------------------"
   echo "Creating conda env $ENV and installing packages"
   echo "----------------------------------------------------------------------"
   conda create -y -n $ENV
   conda activate $ENV
fi

conda install -y -c conda-forge \
              ipython matplotlib meson-python ninja nose notebook numpy \
              petsc4py pip scipy seaborn six

# Build clawpack
echo "----------------------------------------------------------------------"
echo "Building clawpack"
echo "----------------------------------------------------------------------"
pip install --user --no-build-isolation -e .

echo "----------------------------------------------------------------------"
echo "Installation successful, try to run a test"
echo "   conda activate claw"
echo "   cd /tmp"
echo "   python $CLAW/pyclaw/examples/advection_1d/advection_1d.py iplot=1"
echo "It is better not to run inside $CLAW or to modify the code there."
echo "----------------------------------------------------------------------"
