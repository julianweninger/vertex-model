# =============================================================================
# This is a ubuntu-based setup of the Utopia framework.
# It installs all necessary dependencies using ubuntu, clones the repository,
# and sets it up in working condition.
# =============================================================================

ARG BASE_IMAGE="ubuntu:20.04"
FROM ${BASE_IMAGE}

LABEL maintainer="Julian Weninger <julian.weninger@unige.ch>"


# .. Install dependencies .....................................................
# Get ready ...
RUN apt-get update && apt-get install -y hello

# set configurations to default
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y keyboard-configuration
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata

# Install Utopia framework dependencies
RUN apt install -y cmake gcc g++ gfortran git libarmadillo-dev
RUN apt install -y clang clang-tools
RUN apt install -y libboost-all-dev libhdf5-dev libspdlog-dev
RUN apt install -y libyaml-cpp-dev pkg-config
RUN apt install -y python3-dev python3-pip python3-venv


# Optional dependencies
RUN apt-get install -y doxygen
RUN apt-get install -y ffmpeg

# .. Set up Utopia ............................................................
WORKDIR /home
RUN mkdir utopia
WORKDIR /home/utopia

ARG UTOPIA_REPO="https://gitlab.com/utopia-project/utopia.git"
ARG UTOPIA_BRANCH="master"

RUN    git clone ${UTOPIA_REPO} \
    && cd utopia \
    && git checkout ${UTOPIA_BRANCH} \
    && mkdir -p build

WORKDIR /home/utopia/utopia/build

ENV CC=gcc CXX=g++
RUN cmake -DCMAKE_BUILD_TYPE=Release ..

# Can now make the model and perform a test run
RUN pwd
RUN make dummy && ./run-in-utopia-env utopia run dummy


# .. Get ready for using Utopia ...............................................
# TODO Could install the ipython stuff here as well ...

# Communicate to cmake where to find UtopiaConfig.cmake
ENV Utopia_DIR=/home/utopia/utopia/build

# As entrypoint, enter the shell from within the utopia-env
WORKDIR /home
ENTRYPOINT [ "/home/utopia/utopia/build/run-in-utopia-env", \
             "/bin/bash" ]
