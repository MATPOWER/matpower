# Builds a Docker image for Octave 4.2.2 and MATPOWER
#
# Authors:
# Xiangmin Jiao <xmjiao@gmail.com>
# Ray Zimmerman <rz10@cornell.edu>

FROM compdatasci/octave-desktop:latest
LABEL maintainer "Ray Zimmerman <rz10@cornell.edu>"

USER root
WORKDIR /tmp

# Install system packages and Octave
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        coinor-libipopt-dev \
        pkg-config && \
    apt-get clean && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install SDPT3
RUN git clone https://github.com/sqlp/sdpt3.git /opt/sdpt3 && \
    rm -rf /opt/sdpt3/.git && \
    rm -rf /opt/sdpt3/Solvers/Mexfun/*.mex* && \
    rm -rf /opt/sdpt3/Solvers/Mexfun/o_win && \
    octave-cli --no-gui --eval "addpath('/opt/sdpt3', '-end'); savepath" && \
    octave-cli --no-gui --eval "install_sdpt3('-rebuild'); savepath" && \
    rm -rf /opt/sdpt3/Solvers/Mexfun/*.o

# Install SeDuMi
RUN git clone https://github.com/sqlp/sedumi.git /opt/sedumi && \
    rm -rf /opt/sedumi/.git && \
    rm -rf /opt/sedumi/*.mex* && \
    mv /opt/sedumi/vec.m /opt/sedumi/vec.m.disabled && \
    octave-cli --no-gui --eval "addpath('/opt/sedumi', '-end'); savepath" && \
    octave-cli --no-gui --eval "install_sedumi('-rebuild'); savepath" && \
    rm -rf /opt/sedumi/*.o

# Install YALMIP
ENV YALMIP_VER=R20180817
RUN git clone -b ${YALMIP_VER} https://github.com/yalmip/YALMIP.git /opt/yalmip && \
    rm -rf /opt/yalmip/.git && \
    octave-cli --no-gui --eval "addpath(genpath('/opt/yalmip', '.git', 'o_win', 'dev'), '-end'); savepath"

# Install Ipopt
ENV IPOPT_VER=3.11.9
COPY ./.travis/Makefile /opt
# RUN curl -SL https://raw.githubusercontent.com/MATPOWER/matpower/master/.travis/Makefile -o /opt/Makefile
RUN curl -SL https://github.com/coin-or/Ipopt/archive/releases/${IPOPT_VER}.tar.gz | tar -xzC /opt && \
    mv /opt/Ipopt-releases-${IPOPT_VER}/Ipopt/contrib/MatlabInterface /opt/ipopt && \
    mv /opt/Makefile /opt/ipopt/src && \
    make -C /opt/ipopt/src && \
    rm -rf /opt/Ipopt-${IPOPT_VER} && \
    mv /opt/ipopt/src/*.mex /opt/ipopt/ && \
    rm -rf /opt/ipopt/src/*.o && \
    octave-cli --no-gui --eval "addpath('/opt/ipopt', '-end'); savepath"

# Install OSQP
RUN export OCTAVE_VER=`octave-cli --no-gui --eval "fprintf('%s', ver('octave').Version)"` && \
    git clone --recursive https://github.com/oxfordcontrol/osqp /opt/osqp && \
    git clone --recurse-submodules https://github.com/oxfordcontrol/osqp-matlab /opt/osqp-matlab && \
    mkdir /opt/osqp/build && \
    cd /opt/osqp/build && \
    cmake -G "Unix Makefiles" .. && \
    cmake --build . && \
    ln -s /usr/include/octave-$OCTAVE_VER/octave/Matrix.h /opt/matrix.h && \
    ls -al /opt && \
    cd /opt/osqp-matlab && \
    octave-cli --no-gui --eval "mex -I/opt  -I/opt/osqp/include -I/opt/osqp/lin_sys/direct/qdldl -I/opt/osqp/lin_sys/direct/qdldl/qdldl_sources/include osqp_mex.cpp /opt/osqp/build/out/libosqp.a" && \
    ls -al && \
    octave-cli --no-gui --eval "addpath('/opt/osqp-matlab', '-end'); savepath" && \
    rm -rf /opt/osqp /opt/osqp-matlab/.git
COPY ./mp-opt-model/.github/osqp/*.m /opt/osqp-matlab/

# Install test script
COPY ./docker/matpower_docker_tests.sh /usr/local/bin

# Install MATPOWER
COPY . /opt/matpower
RUN git clone https://github.com/MATPOWER/matpower-extras.git /opt/matpower/extras && \
    rm -rf /opt/matpower/extras/.git && \
    octave-cli --no-gui -p /opt/matpower --eval "install_matpower(1,1,1)"

# RUN git clone --depth=5 https://github.com/MATPOWER/matpower.git /opt/matpower && \
#     git clone --depth=5 https://github.com/MATPOWER/matpower-extras.git /opt/matpower/extras && \
#     rm -rf /opt/matpower/.git && \
#     rm -rf /opt/matpower/extras/.git && \
#     octave-cli --no-gui -p /opt/matpower --eval "install_matpower(1,1,1)"

RUN cp ~/.octaverc $DOCKER_HOME

WORKDIR $DOCKER_HOME
