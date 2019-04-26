FROM ubuntu:xenial

MAINTAINER Jaysen <jaysen.sawmynaden@upmc.fr>

RUN apt-get update && \
    apt-get install -qy \
       python \
       csh \
       flex \
       patch \
       gfortran \
       g++ \
       make \
       xorg-dev \
       bison \
       libbz2-dev \
       openmpi-bin \
       cmake=3.5.1-1ubuntu3 \
       libopenmpi-dev



################################################################################

# Set timezone
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Useful tools when running Docker in interactive mode
  # Removing vim-tiny is Ubuntu-specific:
RUN apt-get remove -qy vim-tiny && \
    apt-get -qy install \
      vim \
      htop \
      screen \
      bash-completion \
      wget

#install cmake
#RUN  wget http://www.cmake.org/files/v3.0/cmake-3.0.2.tar.gz
#RUN tar xzf cmake-3.0.2.tar.gz
#RUN cd cmake-3.0.2
#RUN ./configure --prefix=/opt/cmake
#RUN make
#RUN make install

# Extend bash history
RUN sed -i 's/HISTSIZE\=1000/HISTSIZE\=1000000/' /root/.bashrc && sed -i 's/HISTFILESIZE\=2000/HISTFILESIZE\=2000000/' /root/.bashrc
# Modify .bashrc for (improved) autocompletion
RUN sed -i '/^#.*bash_completion/s/^#//' /root/.bashrc && sed -i '$ s/^#//' /root/.bashrc
# Change the default shell in screen to bash
RUN echo "shell \"/bin/bash\"" > /root/.screenrc

# Vim: default syntax highlighting + highlight search
RUN echo "colorscheme default" > /root/.vimrc
RUN echo "set hlsearch" >> /root/.vimrc

################################################################################

# Install GROMACS (version 5.1.2)

COPY scripts/install_gromacs.sh ./

RUN chmod 744 install_gromacs.sh
RUN mkdir -p /home/REMD/data /home/REMD/output/
RUN mkdir -p /home/REMD/src/lunch_REMD/
RUN mkdir -p /home/REMD/src/acpype/
RUN mkdir -p /home/REMD/src/scripts/
RUN mkdir -p /home/REMD/src/scwrl3_lin



COPY ./scripts/lunch_REMD/*.py /home/REMD/src/lunch_REMD/
COPY ./parameters/*.mdp /home/REMD/src/
COPY ./src/acpype/* /home/REMD/src/acpype/
#RUN ./install_gromacs.sh

################################################################################

# Install AmberTools (version 18)

  # Download from RPBS OwnCloud
RUN wget https://owncloud.rpbs.univ-paris-diderot.fr:443/owncloud/index.php/s/5yoyGkC9bbadNJ0/download && mv download amber18.tar.gz

#RUN tar -xzfv amber18.tar.gz
COPY scripts/install_amber.sh ./
RUN chmod 744 install_amber.sh
#RUN ./install_amber.sh

#RUN rm amber18.tar.gz

################################################################################

# Cleanup
RUN apt-get autoremove -qy && \
apt-get clean && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
