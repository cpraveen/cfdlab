VER=$1

sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-$VER $VER
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-$VER $VER
sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-$VER $VER

sudo update-alternatives --config gcc && \
sudo update-alternatives --config g++ && \
sudo update-alternatives --config gfortran
