apt-get upgrade -y
apt-get install build-essential cmake wget openssh-server -y
wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-5.1.2.tar.gz
tar zxvf gromacs-5.1.2.tar.gz 
cd gromacs-5.1.2
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON
make 
make install
cd
source /usr/local/gromacs/bin/GMXRC

