# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:2.3.1-ubuntu24.04
USER root

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    wget \
    nano \
    git \
    patchelf \
    #tcsh\
    #ninja-build \
    #meson \
    libxrender1\
    libxext6 \
    g++ \
    gcc \
    libc6-dev \
    make \
    postgresql \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*  



COPY ./app /home/screener
#RUN chmod +x /home/screener/install/build_deps.sh
RUN chmod +x /home/screener/screener_dock.sh



WORKDIR /home/screener

#create env vars
RUN export APP_PATH=/home/screener
RUN export FILE_PATH=/home/screener/files
RUN export MODULE_PATH=/home/screener/modules
RUN export OPENBLAS_NUM_THREADS=1
RUN export OMP_NUM_THREADS=1
RUN export OPENMM_CPU_THREADS=1
RUN export MKL_NUM_THREADS=1


#save env vars
#RUN echo "$APP_PATH" >> /home/screener/env/.env
#RUN echo "$FILE_PATH" >> /home/screener/env/.env
#RUN echo "$MODULE_PATH" >> /home/screener/env/.env


#create conda env from screener


RUN micromamba create -f /home/screener/env/screener.yml
#RUN micromamba create -f /home/screener/modules/pxdock/pxdock.yml

#RUN micromamba run -n pxdock pip install --use-pep517 git+https://github.com/bytedance/Protenix-Dock.git

#patch executable to work in docker
#RUN patchelf --set-rpath '$ORIGIN:/opt/conda/envs/pxdock/lib' /opt/conda/envs/pxdock/lib/python3.10/site-packages/pxdock/_protenix_dock*.so
#RUN patchelf --set-rpath '$ORIGIN:/opt/conda/envs/pxdock/lib' /opt/conda/envs/pxdock/lib/python3.10/site-packages/pxdock/libdock.so


RUN wget -O /home/screener/install/files/crest.tar.xz --no-check-certificate  https://github.com/crest-lab/crest/releases/download/v3.0.2/crest-gnu-12-ubuntu-latest.tar.xz
RUN tar -xf /home/screener/install/files/crest.tar.xz -C /home/screener/bin
RUN chmod +x /home/screener/bin/crest/crest


RUN micromamba create -f /home/screener/modules/CCTBX/cctbx.yml
RUN wget -P /home/screener/install/files/monomer --no-check-certificate --content-disposition https://github.com/cctbx/cctbx_project/releases/download/v2025.10/chem_data-2025.10-pyhe0d8492_0.conda
RUN micromamba install -n CCTBX /home/screener/install/files/monomer/chem_data-2025.10-pyhe0d8492_0.conda


#RUN micromamba env -n Screener export > /home/screener/app/env/screener.yml

#install protenix




#RUN micromamba create -f /home/screener/modules/protenix/protenix.yml
#RUN micromamba run -n protenix pip install git+https://github.com/bytedance/Protenix.git





#build missing dependencies
#RUN /home/screener/install/build_deps.sh

#flush install files
RUN rm -rf /home/screener/install/files

#get protenix dependencies
#RUN git clone -b v3.5.1 https://github.com/NVIDIA/cutlass.git /opt/cutlass
#ENV CUTLASS_PATH=/opt/cutlass

CMD ["/bin/bash"]       
#CMD  micromamba run -n Screener  python ./scripts/screener.py






