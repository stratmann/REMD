FROM ubuntu
MAINTAINER Jaysen <jaysen.sawmynaden@upmc.fr>
RUN apt-get upgrade -y
RUN apt-get update && apt-get install -y python && apt-get install -y vim
RUN apt-get install -y csh flex patch gfortran g++ make xorg-dev bison libbz2-dev
RUN apt-get install -y openmpi-bin libopenmpi-dev
RUN apt-get install -y wget
COPY install_gromacs.sh ./
COPY ./amber18.tar ./
COPY install_amber.sh ./
RUN chmod 744 install_gromacs.sh install_amber.sh
RUN mkdir -p /home/REMD
RUN mkdir -p /home/REMD/src /home/REMD/data /home/REMD/output
COPY lunch_REMD.py /home/REMD/src/
COPY *.mdp /home/REMD/src/
COPY ./acpype /home/REMD/src/
RUN tar xfv amber18.tar
RUN rm amber18.tar
#RUN ./install_gromacs.sh
#RUN ./install_amber.sh

