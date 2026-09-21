# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:2.3.1-ubuntu24.04

USER root
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    wget \
    nano \
    git 

COPY ./app /home/crest
COPY ./run_simulation.sh /home/crest/run_simulation.sh
RUN chmod +x /home/crest/run_simulation.sh

USER $MAMBA_USER


#create env vars
ENV APP_PATH=/home/crest
ENV SHARE_FILE_PATH=/home/crest/shared_files
ENV SCRIPT_PATH=/home/crest/scripts
ENV PATH="/usr/local/gromacs/bin:${PATH}"
ENV OMP_NUM_THREADS=1,1

RUN micromamba create -f /home/crest/env/crest.yml

WORKDIR /home/crest
CMD ["/bin/bash"]     