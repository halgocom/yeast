# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:2.3.1-ubuntu24.04
USER root

ARG GMX_VERSION=2026.3

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    wget \
    nano \
    git \
    gcc\ 
    cmake\
    libopenmpi-dev\ 
    openmpi-bin


COPY ./app /home/gromacs
COPY ./run_simulation.sh /home/gromacs/run_simulation.sh
RUN chmod +x /home/gromacs/run_simulation.sh

USER $MAMBA_USER
#Softaware Compilation
WORKDIR /home/gromacs/app/install

#GROMACS
#Downaload and unpack 
RUN wget -O ./gromacs.tar.gz "https://ftp.gromacs.org/gromacs/gromacs-$GMX_VERSION.tar.gz"
RUN mkdir ./gromacs_src
RUN tar xfz ./gromacs.tar.gz -C ./gromacs_src
RUN mv ./gromacs_src/gromacs-$GMX_VERSION ./gromacs_src/gromacs

#Install (Takes some time)
WORKDIR /home/gromacs/app/install/gromacs_src/gromacs
RUN mkdir ./build 
WORKDIR /home/gromacs/app/install/gromacs_src/gromacs/build
RUN cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_MPI=OFF -DGMX_OPENMP=ON -DGMXAPI=ON -DBUILD_SHARED_LIBS=ON
RUN make
RUN make check

USER root
RUN make install
RUN source /usr/local/gromacs/bin/GMXRC
USER $MAMBA_USER
#flush installation dir
RUN rm -rf /home/gromacs/app/install/gromacs_src 

#create env vars
ENV APP_PATH=/home/gromacs
ENV SHARE_FILE_PATH=/home/shared_files
ENV SCRIPT_PATH=/home/gromacs/scripts
ENV PATH="/usr/local/gromacs/bin:${PATH}"
ENV gmxapi_ROOT=/usr/local/gromacs

RUN source /usr/local/gromacs/bin/GMXRC

RUN micromamba create -n conda_lock  conda-lock
RUN micromamba create -n openmm python=3.12 -c conda-forge openmm mdtraj matplotlib scipy prody rdkit

WORKDIR /home/gromacs
CMD ["/bin/bash"]     