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
       libopenmpi-dev \
       python-pip
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

#RUN apt-get upgrade -y
RUN apt-get install build-essential cmake wget openssh-server -y
RUN wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-5.1.2.tar.gz
RUN tar zxvf gromacs-5.1.2.tar.gz -C ./
WORKDIR gromacs-5.1.2
RUN mkdir ./build
WORKDIR build
RUN cmake .. -DGMX_BUILD_OWN_FFTW=ON
RUN make
RUN make install
WORKDIR ..
RUN mkdir -p /home/MD/data /home/MD/parameters/FF/ /home/MD/Example/ /home/MD/src


COPY ./scripts/* /home/MD/scripts/
COPY ./parameters/*.mdp /home/MD/parameters/
COPY ./parameters/FF /gromacs-5.1.2/share/top/
COPY ./Example/* /home/MD/Example/
COPY ./src/scwrl3_lin.tar.gz /home/MD/src/
RUN tar -zxvf /home/MD/src/scwrl3_lin.tar.gz -C /home/MD/src/
WORKDIR /home/MD/src/scwrl3_lin/
RUN ./setup
WORKDIR /

WORKDIR /home/MD
################################################################################
# Cleanup
RUN apt-get autoremove -qy && \
apt-get clean && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
