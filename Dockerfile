FROM debian:stable as builder

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -yqq --no-install-recommends \
	ca-certificates \
	less \
	git \
	make \
	pkg-config \
	g++ \
	gcc \
	libc6-dev \
	octave \
	octave-optim \
	liboctave-dev \
	libopenblas-dev \
	coinor-libipopt-dev \
        libglpk-dev \
	&& rm -rf /var/lib/apt/lists/*

ENV CFLAGS -I/usr/include/openblas

RUN git clone https://github.com/yalmip/YALMIP.git /usr/src/yalmip 
RUN git clone https://github.com/sqlp/sdpt3.git /usr/src/sdpt3
RUN git clone https://github.com/sqlp/sedumi.git /usr/src/sedumi
RUN git clone https://github.com/rwl/ipopt.git /usr/src/ipopt
RUN git clone https://github.com/MATPOWER/matpower-extras.git /usr/src/extras

ADD . /usr/src/matpower

RUN mv /usr/src/extras /usr/src/matpower/extras


# GNU Octave 4.0.2 does not include ismembc function.
RUN printf 'function varargout = ismembc(varargin)\n  varargout = cell(nargout, 1);\n  [varargout{:}] = ismember(varargin{:});\nendfunction' >> /usr/src/matpower/extras/sdp_pf/ismembc.m


RUN octave-cli --no-gui --eval "addpath(genpath('/usr/src/yalmip', '.git', 'dev'));savepath"

RUN cd /usr/src/sdpt3 && octave-cli --no-gui --eval "install_sdpt3('-rebuild')"
RUN octave-cli --no-gui --eval "addpath(genpath('/usr/src/sdpt3', '.git', 'o_win'));savepath"

RUN cd /usr/src/sedumi && octave-cli --no-gui --eval "install_sedumi('-rebuild')"
RUN octave-cli --no-gui --eval "addpath(genpath('/usr/src/sedumi', '.git', 'o_win'), '-end');savepath"

RUN make -C /usr/src/ipopt
RUN octave-cli --no-gui --eval "addpath(genpath('/usr/src/ipopt', '.git'));savepath"

RUN cd /usr/src/matpower && octave-cli --no-gui --eval "install_matpower(1,1,1)"

WORKDIR /workspace

CMD ["octave-cli", "-q", "-w", "--no-gui"]


FROM debian:stable

RUN apt-get update && apt-get install -yqq --no-install-recommends \
	octave \
	octave-optim \
	libopenblas-base \
	coinor-libipopt1v5 \
	libglpk40 \
	less \
	&& rm -rf /var/lib/apt/lists/*

ENV INSTALL_PATH /usr/local/matpower

COPY --from=builder /usr/src/matpower $INSTALL_PATH
COPY --from=builder /usr/src/yalmip $INSTALL_PATH/yalmip
COPY --from=builder /usr/src/sdpt3 $INSTALL_PATH/sdpt3
COPY --from=builder /usr/src/sedumi $INSTALL_PATH/sedumi
COPY --from=builder /usr/src/ipopt $INSTALL_PATH/ipopt

RUN octave-cli --no-gui --eval "addpath(genpath('$INSTALL_PATH', '.git', 'o_win', 'dev'), '-end');savepath"

WORKDIR /workspace

ENTRYPOINT ["octave-cli", "-q", "-w", "--no-gui"]

