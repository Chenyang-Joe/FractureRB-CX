# =============================================================
# Dockerfile for legacy C++ project FractureRB
# - Based on Debian Buster
# - Installs Bullet 2.83 and OpenVDB 2.2 using make
# =============================================================

FROM debian:buster

# -------------------------------------
# Install system dependencies
# -------------------------------------
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    g++ \
    libtbb-dev \
    zlib1g-dev \
    libeigen3-dev \
    libboost-all-dev \
    libilmbase-dev \
    libopenexr-dev \
    libblosc-dev \
    wget \
    gdb \
    && rm -rf /var/lib/apt/lists/*

# -------------------------------------
# Build Bullet 2.83 from source (with extras)
# -------------------------------------
WORKDIR /tmp
RUN git clone https://github.com/bulletphysics/bullet3.git && \
    cd bullet3 && \
    git checkout 2.83 && \
    mkdir build && cd build && \
    cmake .. \
      -DBUILD_SHARED_LIBS=ON \
      -DBUILD_EXTRAS=ON \
      -DBUILD_BULLET2_DEMOS=OFF \
      -DBUILD_OPENGL3_DEMOS=OFF \
      -DCMAKE_INSTALL_PREFIX=/usr/local && \
    make -j"$(nproc)" && \
    make install && \
    cp Extras/Serialize/BulletWorldImporter/libBulletWorldImporter.so /usr/local/lib/ && \
    cp Extras/Serialize/BulletFileLoader/libBulletFileLoader.so /usr/local/lib/ && \
    cp -r ../Extras/Serialize/BulletWorldImporter/ /usr/local/include && \
    cp -r ../Extras/Serialize/BulletFileLoader/ /usr/local/include

# -------------------------------------
# Build OpenVDB 2.2 using Makefile (only core lib)
# -------------------------------------
WORKDIR /tmp

RUN wget https://github.com/oneapi-src/oneTBB/archive/refs/tags/2017_U7.tar.gz -O tbb-2017_U7.tar.gz && \
    tar -xzf tbb-2017_U7.tar.gz && \
    cd oneTBB-2017_U7 && \
    make -j$(nproc) 


WORKDIR /tmp
RUN git clone https://github.com/AcademySoftwareFoundation/openvdb.git && \
    cd openvdb && \
    git checkout f9eca01833739a6c24b2fe12a8c3a22ed38e66cc
    # cd openvdb && 
    # make -j"$(nproc)" && \
    # make install
    # docker cp /Users/chenyangxu/Codebase/AI4CE/assembly2025spring/FractureRB-CX/dockertools/Makefile_no_python_no_viewer cf04bb981b31:/tmp/openvdb/openvdb/Makefile \

WORKDIR /tmp/openvdb/openvdb
COPY dockertools/Makefile_no_python_no_viewer Makefile

RUN apt-get update && apt-get install -y libcppunit-dev && \
    apt-get update && apt-get install -y libjemalloc-dev && \
    make HALF_LIB_DIR=/usr/lib/x86_64-linux-gnu \
        TBB_LIB_DIR=/tmp/oneTBB-2017_U7/build/linux_intel64_gcc_cc*_release \
        TBB_INCL_DIR=/tmp/oneTBB-2017_U7/include \
        EXR_LIB_DIR=/usr/lib/x86_64-linux-gnu \
        -j$(nproc) &&\
    apt-get update && apt-get install -y doxygen &&\
    make install INSTALL_DIR=/usr/local

#     sed -i '/python/d' Makefile && \ 
#     rm -rf python unittest viewer && \
#     make libopenvdb.so -j"$(nproc)" && \
#     make install

# -------------------------------------
# fftw, hyena might use it
# -------------------------------------
WORKDIR /tmp
RUN wget https://www.fftw.org/fftw-3.3.10.tar.gz && \
    tar -xzf fftw-3.3.10.tar.gz

WORKDIR /tmp/fftw-3.3.10
RUN apt-get update && apt-get install -y m4 autoconf libtool bash && \
    chmod +x configure && \
    bash ./configure --enable-shared --prefix=/usr/local && \
    make -j4 && \
    make install




# -------------------------------------
# Set up application workspace
# -------------------------------------
WORKDIR /app
COPY . .

RUN apt-get update && apt-get install -y libtclap-dev
RUN apt update && apt install bc
# RUN apt-get install gdb

RUN mkdir -p build
WORKDIR /app/build

# RUN apt-get install gdb
# RUN cmake .. \
#         -DCMAKE_CXX_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=0" \
#         -DCMAKE_BUILD_TYPE=Debug \
#         -DHLIB_INC=/usr/include/eigen3 \
#         -DBULLET_IMPORTER=/usr/local/lib/libBulletWorldImporter.so \
#         -DBULLET_LOADER=/usr/local/lib/libBulletFileLoader.so \
#         -DCMAKE_CXX_FLAGS="-I/tmp/bullet3/Extras/Serialize" \
#         -DHYENA_LIB2=/app/hyena/libHyENAlib2.so \
#         -Dzlib=/usr/lib/x86_64-linux-gnu/libz.so \
#         -DHalflib=/usr/lib/x86_64-linux-gnu/libHalf.so \
#         -DOpenVDBlib=/tmp/openvdb/openvdb/libopenvdb.so \
#         -DOpenVDBinclude=/tmp/openvdb 

# RUN make -j$(nproc)

RUN ls /app/build


# -------------------------------------
# Default: open shell for manual build
# -------------------------------------
CMD ["bash"]
