# Builds a Docker image for Octave and MATPOWER (7+)
# Use the following --build-arg values:
#   BASE_IMAGE=<image> (default is matpower/octave)
#   BASE_TAG=<tag> e.g. for a specific version of Octave (default is latest)
#   MP_SRC=local - to copy the local files to install MATPOWER, pulls
#                  matpower-extras from GitHub
#          github - (default) to clone matpower and matpower-extras from GitHub
#   BRANCH=<branch or tag> to specify the branch or tag for GitHub checkouts
#
# docker build -f docker/Dockerfile-matpower -t matpower/matpower:latest .
# docker build --build-arg BRANCH=7.1 --build-arg BASE_TAG=5.2.0 -f docker/Dockerfile-matpower -t matpower/matpower:7.1-oct-5.2.0 .
#
# Authors:
# Xiangmin Jiao <xmjiao@gmail.com>
# Ray Zimmerman <rz10@cornell.edu>

ARG BASE_IMAGE=matpower/octave
ARG BASE_TAG=latest
ARG MP_SRC=local
ARG BRANCH=master
ARG WARN=0              # set to 1 for Octave 6.x+ and MATPOWER 7.0
ARG PATCH_MOST_TEST=0   # set to 1 for Octave 7.x+ and MATPOWER 7.1
FROM ${BASE_IMAGE}:${BASE_TAG} AS base
LABEL maintainer "Ray Zimmerman <rz10@cornell.edu>"
ARG BRANCH
ARG MP_SRC
ARG WARN
ARG PATCH_MOST_TEST

USER root
ENV HOME /root
WORKDIR /tmp

# Install local MATPOWER files
FROM base AS matpower_local
ONBUILD ARG MP_SRC
ONBUILD ARG BRANCH
ONBUILD COPY ./docker/matpower_docker_tests/matpower_docker_tests_7.sh /usr/local/bin/matpower_docker_tests.sh
ONBUILD COPY . $HOME/packages/matpower
ONBUILD RUN git clone -b ${BRANCH} --depth=1 https://github.com/MATPOWER/matpower-extras.git $HOME/packages/matpower/extras && \
    rm -rf $HOME/packages/matpower/extras/.git && \
    octave-cli --no-gui -p $HOME/packages/matpower --eval "install_matpower(1,1,1)"

# Install MATPOWER files from GitHub
FROM base AS matpower_github
ONBUILD ARG MP_SRC
ONBUILD ARG BRANCH
ONBUILD COPY ./docker/matpower_docker_tests/matpower_docker_tests_7.sh /usr/local/bin/matpower_docker_tests.sh
ONBUILD RUN git clone -b ${BRANCH} --depth=1 https://github.com/MATPOWER/matpower.git $HOME/packages/matpower && \
    git clone -b ${BRANCH} --depth=1 https://github.com/MATPOWER/matpower-extras.git $HOME/packages/matpower/extras && \
    rm -rf $HOME/packages/matpower/.git && \
    rm -rf $HOME/packages/matpower/extras/.git && \
    octave-cli --no-gui -p $HOME/packages/matpower --eval "install_matpower(1,1,1)"

# Select source for MATPOWER files
FROM matpower_${MP_SRC} AS matpower

# conditionally disable a warning in .octaverc (for Octave 6.x+ and MATPOWER 7.0)
FROM matpower AS warn_0
FROM matpower AS warn_1
ONBUILD RUN echo "warning('off', 'Octave:empty-index');" >> $HOME/.octaverc
FROM warn_${WARN} as warn

FROM warn AS patch_0
FROM warn AS patch_1
ONBUILD COPY ./docker/patches/t_most_mpopf.patch $HOME/packages/
ONBUILD RUN patch -u -b $HOME/packages/matpower/most/lib/t/t_most_mpopf.m -i $HOME/packages/t_most_mpopf.patch
FROM patch_${PATCH_MOST_TEST}

WORKDIR "${HOME}"
