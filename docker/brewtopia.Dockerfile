# =============================================================================
# This is a homebrew-based setup of the Utopia framework.
# It installs all necessary dependencies using homebrew, clones the repository,
# and sets it up in working condition.
# =============================================================================

ARG BASE_IMAGE="homebrew/ubuntu20.04:latest"
FROM ${BASE_IMAGE}

LABEL maintainer="Julian Weninger <julian.weninger@unige.ch>"


# .. Install dependencies .....................................................
# Get ready ...
RUN brew update && brew install hello

# Install Utopia framework dependencies (except boost, see below) and some
# additional dependencies
RUN brew install --display-times gcc llvm pkg-config cmake python@3.9
RUN brew install --display-times armadillo boost hdf5 fmt spdlog yaml-cpp
RUN brew install --display-times ffmpeg graphviz doxygen texlive fftw

# Make sure the correct gcc and Python are linked
RUN brew link gcc
RUN brew link python@3.9 --force

# Use an earlier boost version (1.71) ... from apt :grimacing:
# RUN    apt-get update \
#     && apt-get install -y libboost-dev libboost-test-dev \
#     && apt-get clean
# NOTE There seem to be some issues with boost v1.74, so this is a workaround
#      to get an earlier version installed. (Much easier than via brew.)


# .. Set up Utopia ............................................................
RUN mkdir utopia
WORKDIR /home/linuxbrew/utopia

ARG UTOPIA_CLONE_URL="https://ts-gitlab.iup.uni-heidelberg.de/utopia/utopia.git"
ARG UTOPIA_BRANCH="master"

RUN    git clone ${UTOPIA_CLONE_URL} \
    && cd utopia \
    && git checkout ${UTOPIA_BRANCH} \
    && mkdir -p build

WORKDIR /home/linuxbrew/utopia/utopia/build

ENV CC=gcc CXX=g++
RUN    cmake -DCMAKE_BUILD_TYPE=Release -DHDF5_ROOT=$(brew --prefix hdf5) .. \
    && make dummy


# .. Get ready for using Utopia ...............................................
# TODO Could install the ipython stuff here as well ...

# As entrypoint, enter the shell from within the utopia-env
WORKDIR /home/linuxbrew
ENTRYPOINT [ "/home/linuxbrew/utopia/utopia/build/run-in-utopia-env", \
             "/bin/bash" ]

