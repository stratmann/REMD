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


RUN apt-get install build-essential cmake wget openssh-server -y
RUN wget http://ftp.gromacs.org/pub/gromacs/gromacs-5.1.2.tar.gz
RUN tar zxvf gromacs-5.1.2.tar.gz -C ./
WORKDIR gromacs-5.1.2
RUN mkdir ./build
WORKDIR build
RUN cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_MPI=on
RUN make
RUN make install
WORKDIR /


RUN mkdir -p /home/REMD/data /home/REMD/output/
RUN mkdir -p /home/REMD/scripts/lunch_REMD/ /home/REMD/scripts/analyse_REMD
RUN mkdir -p /home/REMD/src/acpype/
RUN mkdir -p /home/REMD/src/scripts/

COPY ./data/seq*.txt /home/REMD/data/
COPY ./data/RGDpV.pdb /home/REMD/data/RGDpV.pdb

COPY ./scripts/lunch_REMD/*.py /home/REMD/scripts/lunch_REMD/
COPY ./scripts/analyse_REMD/*.py /home/REMD/scripts/analyse_REMD/
COPY ./parameters/*.mdp /home/REMD/src/
COPY ./src/acpype/* /home/REMD/src/acpype/
RUN chmod 744 /home/REMD/src/acpype/acpype.py
COPY ./src/scwrl3_lin.tar.gz /home/REMD/src/
RUN tar -zxvf /home/REMD/src/scwrl3_lin.tar.gz -C /home/REMD/src/
COPY ./src/BackboneReference /home/REMD/src/scwrl3_lin
WORKDIR /home/REMD/src/scwrl3_lin/
#RUN ./setup
WORKDIR /




################################################################################

# Install AmberTools (version 18)

  # Download from RPBS OwnCloud
RUN wget https://owncloud.rpbs.univ-paris-diderot.fr:443/owncloud/index.php/s/5yoyGkC9bbadNJ0/download && mv download amber18.tar.gz

RUN tar -zxvf amber18.tar.gz -C /
WORKDIR /amber18
RUN export AMBERHOME=`pwd`
#######A automatiser ces lignes de commandes########
RUN yes | ./configure gnu
RUN /bin/bash -c "source /amber18/amber.sh"
RUN make install
RUN echo "source /amber18/amber.sh" >> ~/.bashrc 

#RUN chmod 744 install_amber.sh

RUN rm /amber18.tar.gz
WORKDIR /home/REMD/

################################################################################



#Une fois que ambertools a été installé
#RUN /amber18/miniconda/bin/conda install -c omnia mdtraj pyemma -y
#RUN python -mpip install -U pip
#RUN python -mpip install -U matplotlib
#RUN python -mpip install numpy mdtraj 
#RUN python -mpip install pyemma


################################################################################

# Cleanup
RUN apt-get autoremove -qy && \
apt-get clean && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

#Execute 
WORKDIR /home/REMD/src/scwrl3_lin/
RUN ./setup
WORKDIR /home/REMD/
#sudo docker pull continuumio/anaconda
#sudo docker run -it --rm continuumio/anaconda /bin/bash
#sudo conda install -c omnia pyemma

# Add sudo
#RUN apt-get -y install sudo

# Add user ubuntu with no password, add to sudo group
#RUN adduser --disabled-password --gecos '' ubuntu
#RUN adduser ubuntu sudo
#RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
#USER ubuntu
#WORKDIR /home/ubuntu/
#RUN chmod a+rwx /home/ubuntu/
#RUN echo `pwd`

# Anaconda installing
RUN wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
RUN bash Anaconda3-5.0.1-Linux-x86_64.sh -b
RUN rm Anaconda3-5.0.1-Linux-x86_64.sh

# Set path to conda
ENV PATH /root/anaconda3/bin:$PATH
ENV PATH /home/ubuntu/anaconda3/bin:$PATH
RUN echo "y" | conda install -c omnia pyemma
################################################################################
####Pour régler les soucis avec matplotlib
RUN apt-get update
RUN apt-get install -y libgl1-mesa-dev
