FROM ubuntu:20.04 
LABEL maintainer="borjafreire1" 
LABEL version="1.0" 
LABEL description="ViQUF docker image"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y python3 cmake pip software-properties-common wget nano curl build-essential procps file git && apt-get clean
RUN pip install install numpy sklearn scipy biopython pandas matplotlib plotly altair gurobipy

# GH
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-key C99B11DEB97541F0 \
    && apt-add-repository https://cli.github.com/packages \
    && apt update \
    && apt install gh

# ENV FOLDER_PROJECT /var/borjafreire1
# RUN mkdir -p $FOLDER_PROJECT
# run cd $FOLDER_PROJECT
# COPY install_gh.sh $FOLDER_PROJECT
# run ./$FOLDER_PROJECT/install_gh.sh

# Repositories: ViQUF, Gatb
run git clone https://github.com/borjaf696/ViQUF
run cd ViQUF \
    && cd lib/ \ 
    && rm -r gatb-core \ 
    && git clone https://github.com/GATB/gatb-core \
    && cd gatb-core/gatb-core \
    && mkdir build ; cd build ; cmake .. ; make -j8 
# Lemon
run cd ViQUF/lib/ \
    && wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz \
    && tar xvf lemon-1.3.1.tar.gz \
    && cd lemon-1.3.1 \
    && cd build \
    && cmake .. \
    && make install

# SDSL
 run git clone https://github.com/simongog/sdsl-lite.git \
    && cd sdsl-lite \
    && ./install.sh
    
# Make ViqUF
run cd ViQUF\
    && make clean && make