# docker build -t hosnwt .
# docker run -it hosnwt /bin/bash

FROM debian
LABEL maintainer "guillaume.jacquenot@gmail.com"

RUN apt-get update && \
    apt-get install -y \
        gfortran \
        cmake \
        liblapack-dev \
        fftw3 \
        libfftw3-dev

WORKDIR .
ADD . /hos-nwt
RUN cd /hos-nwt && \
    cd cmake && \
    mkdir -p build && \
    cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make && \
    make test

# RUN cd /hos-nwt && \
#     cd cmake && \
#     mkdir -p build && \
#     cd build && \
#     cmake .. -DCMAKE_BUILD_TYPE=Coverage && \
#     make && \
#     make test && \
#     make coverage
