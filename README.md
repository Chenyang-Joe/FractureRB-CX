# FractureRB

Sources for the paper "Fast approximations for boundary element based brittle fracture simulation".
http://pub.ist.ac.at/group_wojtan/projects/2016_Hahn_FastFracture/



This project replaces https://github.com/david-hahn/FractureBEM.

Your OpenVDB version should include the bugfix to TreeIterator.h from
https://github.com/dreamworksanimation/openvdb/commit/e6aa3c83b58a90520c9825a2b385b8185b5c0df3

---

# My deployment

I deploy this project using Docker on Debian (debian::burst) while developing locally on macOS. The g++ version is 8.3.0-6. The hyena package I use is HyENAlib2_x64_Debian_g++-4.9.zip. According to the debug info (see the last part), the bugs are possibly caused by the hyena compiling environment latency.

### Running the Program

#### How to Run the `glass.sh` Example1. build docker images
1. **Build the Docker Image**
```bash
docker build -t fracture-rb . 
```

2. **Enter the Docker Container**
```bash
docker run -it -v "$PWD":/app \
-w /app \
fracture-rb bash
```


<!-- 3. **Install GDB for Debugging**
```bash
apt-get install gdb
``` -->

3. **Build and Compile the Project**
```bash
mkdir -p build && cd build
cmake .. \
-DCMAKE_CXX_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=0" \
-DCMAKE_BUILD_TYPE=Debug \
-DHLIB_INC=/usr/include/eigen3 \
-DBULLET_IMPORTER=/usr/local/lib/libBulletWorldImporter.so \
-DBULLET_LOADER=/usr/local/lib/libBulletFileLoader.so \
-DCMAKE_CXX_FLAGS="-I/opt/bullet3/Extras/Serialize" \
-DHYENA_LIB2=/app/hyena/libHyENAlib2.so \
-Dzlib=/usr/lib/x86_64-linux-gnu/libz.so \
-DHalflib=/usr/lib/x86_64-linux-gnu/libHalf.so \
-DOpenVDBlib=/usr/local/lib/libopenvdb.so \
-DOpenVDBinclude=/usr/local/include

make -j$(nproc)
```

4. **Run the Example**
```bash
cd ../examples
chmod +x glass.sh
./glass.sh
run
bt
```

---

## Debugging Information

With the new libHyENAlib2.so, the glass example starts to generate cracks instead of just crashes directly. However, it still crashes in around 30 mins and I have attached the debug information which is related to HyENA.
```bash
(gdb) bt
#0  0x00007ffff763a0c6 in void hyena::InnerAssembler<hyena::DofTraits<(hyena::ELEMENT_SHAPE)2, (hyena::APPROXIMATION)2, (hyena::APPROXIMATION)2, (hyena::PROBLEM)2, (hyena::SPACE_TYPE)1> >::operator()<__gnu_cxx::__normal_iterator<hyena::Mat<double, 3, 1> const*, std::vector<hyena::Mat<double, 3, 1>, std::allocator<hyena::Mat<double, 3, 1> > > >, Eigen::Matrix<double, -1, 1, 0, -1, 1>, hyena::InnerIntegrator3d<hyena::ColloKernel<hyena::DofTraits<(hyena::ELEMENT_SHAPE)2, (hyena::APPROXIMATION)2, (hyena::APPROXIMATION)2, (hyena::PROBLEM)2, (hyena::SPACE_TYPE)1>, hyena::DofTraits<(hyena::ELEMENT_SHAPE)2, (hyena::APPROXIMATION)2, (hyena::APPROXIMATION)2, (hyena::PROBLEM)2, (hyena::SPACE_TYPE)1>, hyena::Elastostatic3D, hyena::DLP_tag> > >(__gnu_cxx::__normal_iterator<hyena::Mat<double, 3, 1> const*, std::vector<hyena::Mat<double, 3, 1>, std::allocator<hyena::Mat<double, 3, 1> > > >, __gnu_cxx::__normal_iterator<hyena::Mat<double, 3, 1> const*, std::vector<hyena::Mat<double, 3, 1>, std::allocator<hyena::Mat<double, 3, 1> > > >, std::vector<hyena::LDof<(hyena::ELEMENT_SHAPE)2, (hyena::APPROXIMATION)2, (hyena::APPROXIMATION)2, (hyena::PROBLEM)2, (hyena::SPACE_TYPE)1> const*, std::allocator<hyena::LDof<(hyena::ELEMENT_SHAPE)2, (hyena::APPROXIMATION)2, (hyena::APPROXIMATION)2, (hyena::PROBLEM)2, (hyena::SPACE_TYPE)1> const*> > const&, Eigen::Matrix<double, -1, 1, 0, -1, 1> const&, hyena::InnerIntegrator3d<hyena::ColloKernel<hyena::DofTraits<(hyena::ELEMENT_SHAPE)2, (hyena::APPROXIMATION)2, (hyena::APPROXIMATION)2, (hyena::PROBLEM)2, (hyena::SPACE_TYPE)1>, hyena::DofTraits<(hyena::ELEMENT_SHAPE)2, (hyena::APPROXIMATION)2, (hyena::APPROXIMATION)2, (hyena::PROBLEM)2, (hyena::SPACE_TYPE)1>, hyena::Elastostatic3D, hyena::DLP_tag> > const&, Eigen::Matrix<double, -1, 1, 0, -1, 1>&, double) const () from /app/hyena/libHyENAlib2.so
#1  0x00007ffff762dd79 in void hyena::InnerEvaluation<hyena::DofTraits<(hyena::ELEMENT_SHAPE)2, (hyena::APPROXIMATION)2, (hyena::APPROXIMATION)2, (hyena::PROBLEM)2, (hyena::SPACE_TYPE)1>, hyena::DofTraits<(hyena::ELEMENT_SHAPE)2, (hyena::APPROXIMATION)2, (hyena::APPROXIMATION)2, (hyena::PROBLEM)2, (hyena::SPACE_TYPE)1>, hyena::Elastostatic3D, hyena::DLP_tag>::operator()<Eigen::Matrix<double, -1, 1, 0, -1, 1>, std::vector<hyena::Mat<double, 3, 1>, std::allocator<hyena::Mat<double, 3, 1> > > >(std::vector<hyena::LDof<(hyena::ELEMENT_SHAPE)2, (hyena::APPROXIMATION)2, (hyena::APPROXIMATION)2, (hyena::PROBLEM)2, (hyena::SPACE_TYPE)1> const*, std::allocator<hyena::LDof<(hyena::ELEMENT_SHAPE)2, (hyena::APPROXIMATION)2, (hyena::APPROXIMATION)2, (hyena::PROBLEM)2, (hyena::SPACE_TYPE)1> const*> > const&, Eigen::Matrix<double, -1, 1, 0, -1, 1> const&, std::vector<hyena::Mat<double, 3, 1>, std::allocator<hyena::Mat<double, 3, 1> > > const&, double const&, Eigen::Matrix<double, -1, 1, 0, -1, 1>&) ()
   from /app/hyena/libHyENAlib2.so
#2  0x00007ffff7622550 in FractureSim::HyENAWrapperImpl::computeCrackBaseDisplacements() () from /app/hyena/libHyENAlib2.so
#3  0x0000555555bf333e in FractureSim::FractureBEM::writeMesh (this=0x7fffbc3195a0,
    filename="_out/glass_pCube1_F1_F8_F43_1_147_0", visualQuality=2, visDisplace=true, visCOD=true, visClose=false,
    visOBJ=false) at /app/src/FractureBEM.cpp:520
#4  0x0000555555a90e7a in FractureSim::FractureRB::runFractureSim (this=0x7fffac32d9a0, maxTime=0.0040000000000000001,
    rbTimeCode=147) at /app/src/FractureRB.cpp:236
#5  0x0000555555a70d69 in FractureSim::BulletWrapper::stepSimulation (this=0x7fffffffdd00, outPtr=0x7fffffffe120,
    outFmt=FractureSim::BulletWrapper::OUT_MEL, dt=0.0040000000000000001, impulseThreshold=100000, forceThreshold=10000000)
    at /app/src/BulletWrapper.cpp:447
#6  0x0000555555a5d52b in main (argc=11, argv=0x7fffffffeca8) at /app/src/main.cpp:46
```