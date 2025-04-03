# FractureRB

Sources for the paper "Fast approximations for boundary element based brittle fracture simulation".
http://pub.ist.ac.at/group_wojtan/projects/2016_Hahn_FastFracture/



This project replaces https://github.com/david-hahn/FractureBEM.

Your OpenVDB version should include the bugfix to TreeIterator.h from
https://github.com/dreamworksanimation/openvdb/commit/e6aa3c83b58a90520c9825a2b385b8185b5c0df3


# My deployment

I use docker to deploy this project at debian::burst. My local environment is macos.
To run the program:

# How to run example glass.sh?
1. build docker images
docker build -t fracture-rb .

2. enter bash
docker run -it -v "$PWD":/app \
-w /app \
fracture-rb bash

4. install gdb for debug
apt-get install gdb

5. build and make project
mkdir -p build && cd build
cmake .. \
-DCMAKE_CXX_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=0" \
-DCMAKE_BUILD_TYPE=Debug \
-DHLIB_INC=/usr/include/eigen3 \
-DBULLET_IMPORTER=/usr/local/lib/libBulletWorldImporter.so \
-DBULLET_LOADER=/usr/local/lib/libBulletFileLoader.so \
-DCMAKE_CXX_FLAGS="-I/tmp/bullet3/Extras/Serialize" \
-DHYENA_LIB2=/app/hyena/libHyENAlib2.so \
-Dzlib=/usr/lib/x86_64-linux-gnu/libz.so \
-DHalflib=/usr/lib/x86_64-linux-gnu/libHalf.so \
-DOpenVDBlib=/tmp/openvdb/openvdb/libopenvdb.so \
-DOpenVDBinclude=/tmp/openvdb

make -j$(nproc)

6. run example
cd ../examples
chmod +x glass.sh
./glass.sh
run
bt

# debug info
use bt to trace bugs:

(gdb) run
Starting program: /app/build/FractureRB glass.bullet glass.csv -o _out/glass_ -i 1e5 -f 1e7 -s 2
[Thread debugging using libthread_db enabled]
Using host libthread_db library "/lib/x86_64-linux-gnu/libthread_db.so.1".

% parsing command line: 
% /app/build/FractureRB ...
%   glass.bullet glass.csv -o _out/glass_ -i 1e5 -f 1e7 -s 2
% reading scene glass.bullet with parameters glass.csv
unknown chunk
unknown chunk
unknown chunk
unknown chunk

% object 0: "collisionObject0" 
% ... is sleeping
% ... is static, shape type 28
% ... friction 0.500, restitution 0.000
% object 1: "collisionObject1" --> "nurbsSphere1"
% ... building breakable rigid body:
% ... crack mesh size 0.6 / voxel size 0.1 = resolution ratio 6.0
% ... Young's modulus 3.1e+09, Poisson's ratio 0.327, density 8.9e+03
% ... tensile strength 7.6e+99, tensile toughness 1e+03, compressive factor 3.000
% ... built VDB sphere with radius 2 and voxel size 0.1
% ... building BEM mesh from VDB surface[New Thread 0x7ffff6600700 (LWP 360)]
[New Thread 0x7fffee000700 (LWP 361)]
[New Thread 0x7ffff6000700 (LWP 362)]
[New Thread 0x7ffff5a00700 (LWP 363)]
[New Thread 0x7ffff5400700 (LWP 364)]
[New Thread 0x7ffff4e00700 (LWP 365)]
[New Thread 0x7ffff4800700 (LWP 366)]
[New Thread 0x7fffef800700 (LWP 367)]
[New Thread 0x7fffeec00700 (LWP 370)]
[New Thread 0x7fffef200700 (LWP 369)]
[New Thread 0x7fffefe00700 (LWP 368)]
[New Thread 0x7fffee600700 (LWP 371)]
[New Thread 0x7fffece00700 (LWP 373)]
[New Thread 0x7fffeda00700 (LWP 374)]
[New Thread 0x7fffed400700 (LWP 372)]
[New Thread 0x7fffe7e00700 (LWP 375)]
[New Thread 0x7fffe7400700 (LWP 376)]
[New Thread 0x7fffe6a00700 (LWP 377)]
[New Thread 0x7fffe6000700 (LWP 378)]
[New Thread 0x7fffe5600700 (LWP 379)]
[New Thread 0x7fffe4c00700 (LWP 380)]
[New Thread 0x7fffdbe00700 (LWP 381)]
[New Thread 0x7fffdb400700 (LWP 382)]
[New Thread 0x7fffdaa00700 (LWP 383)]
[New Thread 0x7fffda000700 (LWP 384)]
[New Thread 0x7fffd9600700 (LWP 385)]
[New Thread 0x7fffd8c00700 (LWP 386)]
[New Thread 0x7fffd3e00700 (LWP 387)]
[New Thread 0x7fffd3400700 (LWP 388)]
[New Thread 0x7fffd2a00700 (LWP 389)]

% ... have 60 elements in the BEM mesh
% ... using default material model
% ... timestep is 0.00184, Rayleigh wave speed is 326
% ... built breakable rigid body (0.3196s)
% ... is not sleeping
% ... mass is 2.987e+05
% ... is active, shape type 31
% ... friction 0.300, restitution 0.900
% object 2: "collisionObject2" 
% ... is sleeping
% ... is static, shape type 31
% ... friction 1.000, restitution 0.100
% object 3: "collisionObject3" --> "pCube1"
% ... building breakable rigid body:
% ... crack mesh size 0.4 / voxel size 0.04 = resolution ratio 10.0
% ... Young's modulus 3.1e+09, Poisson's ratio 0.327, density 1e+03
% ... tensile strength 7.6e+03, tensile toughness 1e+03, compressive factor 6.000
% ... created coarse tri-mesh from box shape - use remeshing!
% ... initializing implicit surface
% ... building BEM mesh from VDB surface
% ... have 1500 elements in the BEM mesh
% ... using default material model
% ... timestep is 0.000411, Rayleigh wave speed is 973
% ... built breakable rigid body (10.7357s)
% ... is sleeping
% ... mass is 1.998e+05
% ... is active, shape type 31
% ... friction 0.900, restitution 0.200
% split impulse threshold set to 0.08 (2 * min. voxel size)
% scene import done

% timestepping the rigid body sim ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 41: I 2.354e+05   F 1.866e+07  (nurbsSphere1)
% 41: I 2.354e+05   F 1.905e+07  (pCube1)
% ************************************************************
% fracture simulation for nurbsSphere1_41 ... 2 time steps
% ... initializing BEM solver
Thread 1 "FractureRB" received signal SIGSEGV, Segmentation fault.
__rawmemchr_avx2 () at ../sysdeps/x86_64/multiarch/memchr-avx2.S:65
65      ../sysdeps/x86_64/multiarch/memchr-avx2.S: No such file or directory.

(gdb) bt
#0  __rawmemchr_avx2 () at ../sysdeps/x86_64/multiarch/memchr-avx2.S:65
#1  0x00007ffff71b82c2 in _IO_str_init_static_internal (sf=sf@entry=0x7fffffffc420, 
    ptr=ptr@entry=0x8 <error: Cannot access memory at address 0x8>, size=size@entry=0, pstart=pstart@entry=0x0)
    at strops.c:41
#2  0x00007ffff71ab50d in _IO_vsscanf (string=0x8 <error: Cannot access memory at address 0x8>, 
    format=0x7ffff7648987 "%lu(%lf,%lf,%lf)", args=args@entry=0x7fffffffc550) at iovsscanf.c:40
#3  0x00007ffff71a5cc4 in __sscanf (s=<optimized out>, format=<optimized out>) at sscanf.c:32
#4  0x00007ffff76187dc in FractureSim::HyENAWrapperImpl::parseBoundaryCnd(char const*, hyena::BCData&) ()
   from /app/hyena/libHyENAlib2.so
#5  0x00007ffff761a8ef in FractureSim::HyENAWrapperImpl::buildBoundaryCnd(std::vector<std::string, std::allocator<std::string> > const&) () from /app/hyena/libHyENAlib2.so
#6  0x0000555555bee744 in FractureSim::FractureBEM::initBEM (this=0x555555ee67b0, youngsModulus=3100000000, 
    poissonsRatio=0.32700000000000001, is2D=false, bndCnds=std::vector of length 60, capacity 64 = {...})
    at /app/src/FractureBEM.cpp:32
#7  0x0000555555bee670 in FractureSim::FractureBEM::initBEM (this=0x555555ee67b0, material=..., is2D=false, 
    bndCnds=std::vector of length 60, capacity 64 = {...}) at /app/src/FractureBEM.cpp:22
#8  0x0000555555a90b22 in FractureSim::FractureRB::runFractureSim (this=0x555555eef240, 
    maxTime=0.0040000000000000001, rbTimeCode=41) at /app/src/FractureRB.cpp:201
#9  0x0000555555a70d69 in FractureSim::BulletWrapper::stepSimulation (this=0x7fffffffdd00, outPtr=0x7fffffffe120, 
    outFmt=FractureSim::BulletWrapper::OUT_MEL, dt=0.0040000000000000001, impulseThreshold=100000, 
    forceThreshold=10000000) at /app/src/BulletWrapper.cpp:447
#10 0x0000555555a5d52b in main (argc=11, argv=0x7fffffffeca8) at /app/src/main.cpp:46
