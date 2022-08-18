# =============================================================================
# Expects a version of the aptopia image and upgrades the installed
# dependencies as well as the Utopia installation. Finally, it makes
# the setup compatible for use in the GitLab CI/CD for a models repository.
#
# Suggested base image:  julianweninger/aptopia:latest
#
# NOTE If this build starts to take longer and longer, make sure to run the
#      base image build anew.
# =============================================================================

ARG BASE_IMAGE
FROM ${BASE_IMAGE}

LABEL maintainer="Julian Weninger <julian.weninger@unige.ch>"

# Update and upgrade dependencies ...
RUN apt-get update && apt-get install -y hello
RUN apt-get upgrade -y


# ----------------------------------------------------------------------------

# Update Utopia itself, reconfigure, and build
WORKDIR /home/utopia/utopia

ARG UTOPIA_BRANCH="v1"
RUN git config --global http.sslverify false
RUN    git checkout ${UTOPIA_BRANCH} \
    && git pull

ENV CC=gcc CXX=g++ CXX_FLAGS="-Og"
RUN rm -rf build && mkdir -p build

WORKDIR /home/utopia/utopia/build
ENV CC=gcc CXX=g++
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make dummy && ./run-in-utopia-env utopia run dummy

# Done. Overwrite entry point to make it compatible with GitLab CI/CD
WORKDIR /home
ENTRYPOINT [ ]
