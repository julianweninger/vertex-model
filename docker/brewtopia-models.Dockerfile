# =============================================================================
# Expects a version of the brewtopia image and upgrades the installed
# dependencies as well as the Utopia installation. Finally, it makes
# the setup compatible for use in the GitLab CI/CD for a models repository.
#
# Suggested base image:  blsqr/brewtopia:latest
#
# NOTE If this build starts to take longer and longer, make sure to run the
#      base image build anew.
# =============================================================================

ARG BASE_IMAGE
FROM ${BASE_IMAGE}

LABEL maintainer="Julian Weninger <julian.weninger@unige.ch>"

# Update and upgrade dependencies ...
RUN brew update && brew install hello
RUN brew upgrade


# ----------------------------------------------------------------------------

# Update Utopia itself, reconfigure, and build
WORKDIR /home/linuxbrew/utopia/utopia

ARG UTOPIA_BRANCH="master"
RUN git config --global http.sslverify false
RUN    git checkout ${UTOPIA_BRANCH} \
    && git pull

ENV CC=gcc CXX=g++ CXX_FLAGS="-Og"
RUN    rm -rf build && mkdir -p build \
    && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release -DHDF5_ROOT=$(brew --prefix hdf5) .. \
    && make dummy \
    && ./run-in-utopia-env utopia run dummy

# Done. Overwrite entry point to make it compatible with GitLab CI/CD
WORKDIR /home/linuxbrew
ENTRYPOINT [ ]
