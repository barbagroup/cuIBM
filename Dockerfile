# Dockerfile for cuIBM

FROM nvidia/cuda:8.0-devel-ubuntu16.04
MAINTAINER Olivier Mesnard <mesnardo@gwu.edu>

# install base system
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    make \
    vim \
    git \
    g++ \
    wget \
    ca-certificates

# CUSP
RUN CUSP_VERSION=0.5.1 && \
    TARBALL=v${CUSP_VERSION}.tar.gz && \
    wget https://github.com/cusplibrary/cusplibrary/archive/${TARBALL} && \
    CUSP_DIR=/opt/cusp/${CUSP_VERSION} && \
    mkdir -p ${CUSP_DIR} && \
    tar -xvf ${TARBALL} -C ${CUSP_DIR} --strip-components=1 && \
    rm -f ${TARBALL}

ENV CUSP_DIR=/opt/cusp/0.5.1

# Boost
RUN BOOST_VERSION=1.64.0 && \
    TARBALL=boost_1_64_0.tar.gz && \
    BOOST_DIR=/opt/boost/${BOOST_VERSION} && \
    mkdir -p ${BOOST_DIR} && \
    cd ${BOOST_DIR} && \
    wget https://dl.bintray.com/boostorg/release/${BOOST_VERSION}/source/${TARBALL} && \
    tar -xvf ${TARBALL} -C ${BOOST_DIR} --strip-components=1 && \
    rm -f ${TARBALL}

ENV BOOST_DIR=/opt/boost/1.64.0

# cuIBM
RUN CUIBM_DIR=/opt/cuIBM && \
    mkdir -p ${CUIBM_DIR} && \
    cd /opt && \
    git clone https://github.com/mesnardo/cuIBM.git && \
    cd ${CUIBM_DIR} && \
    BRANCH_NAME=master && \
    git checkout -b ${BRANCH_NAME} origin/${BRANCH_NAME} && \
    make -j2

ENV PATH=/opt/cuIBM/bin:$PATH \
    CUIBM_DIR=/opt/cuIBM

# Miniconda
RUN MINICONDA_SCRIPT=Miniconda3-latest-Linux-x86_64.sh && \
    wget https://repo.continuum.io/miniconda/${MINICONDA_SCRIPT} && \
    bash ${MINICONDA_SCRIPT} -b -p /opt/miniconda && \
    PATH=/opt/miniconda/bin:$PATH && \
    conda update conda && \
    conda install numpy scipy matplotlib pandas && \
    rm -f ${MINICONDA_SCRIPT} && \
    mkdir -p $HOME/.config/matplotlib && \
    echo "backend: agg" > $HOME/.config/matplotlib/matplotlibrc

ENV PATH=/opt/miniconda/bin:$PATH

# Setup snake
RUN cd ${CUIBM_DIR}/external/snake-0.3 && \
    python setup.py install

ENV SNAKE=${CUIBM_DIR}/external/snake-0.3
