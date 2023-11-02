# Base image
FROM ubuntu:20.04

# Metadata
LABEL maintainer="borjafreire1" \
      version="1.0" \
      description="ViQUF docker image"

# Environment settings
ARG DEBIAN_FRONTEND=noninteractive

# Install essential packages
RUN apt update && \
    apt install -y \
        python3 \
        pip \
        cmake \
        software-properties-common \
        wget \
        nano \
        curl \
        build-essential \
        procps \
        file \
        g++ \
        doxygen \
        git && \
    apt-get clean

# Install Python packages
RUN pip3 install numpy scikit-learn scipy biopython pandas matplotlib plotly altair gurobipy

# Install GitHub CLI
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-key 23F3D4EA75716059 && \
    apt-add-repository https://cli.github.com/packages && \
    apt update && \
    apt install gh

# Clone repositories
RUN git clone https://github.com/borjaf696/ViQUF.git

# Install gatb-core
RUN cd ViQUF/lib/ && \
    rm -r gatb-core && \
    git clone https://github.com/GATB/gatb-core && \
    cd gatb-core/gatb-core && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j8

# Install Lemon
RUN cd ViQUF/lib/ && \
    wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz && \
    tar xvf lemon-1.3.1.tar.gz && \
    cd lemon-1.3.1 && \
    cd build && \
    cmake .. && \
    make install

# Install SDSL
RUN git clone https://github.com/simongog/sdsl-lite.git && \
    cd sdsl-lite && \
    ./install.sh

# Build ViqUF
RUN cd ViQUF && \
    make clean && \
    make
