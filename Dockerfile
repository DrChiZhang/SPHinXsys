FROM ubuntu:20.04

ARG build_with_dependencies_source=0
ARG sph_only_static_build=0

ENV TZ=Europe/Berlin
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && apt-get install -y \ 
    apt-utils \
    build-essential \
    cmake \
    googletest \
    libtbb-dev \
    libboost-all-dev \
    libsimbody-dev \
    libsimbody3.6 \      
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV TBB_HOME=/usr/lib/x86_64-linux-gnu
ENV BOOST_HOME=/usr/lib/x86_64-linux-gnu
ENV SIMBODY_HOME=/usr

COPY ./ /home/SPHinXsys/
WORKDIR /home/SPHinXsys
RUN rm -rf build
RUN mkdir build && cd build && cmake .. -DBUILD_WITH_DEPENDENCIES_SOURCE=${build_with_dependencies_source} -DSPH_ONLY_STATIC_BUILD=${sph_only_static_build} && make -j$(nproc)