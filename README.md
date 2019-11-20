# High-Performance Material Point Method (CB-Geo mpm)
> [CB-Geo Computational Geomechanics Research Group](https://www.cb-geo.com)

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/cb-geo/mpm/develop/license.md)
[![Developer docs](https://img.shields.io/badge/developer-docs-blue.svg)](http://cb-geo.github.io/mpm)
[![User docs](https://img.shields.io/badge/user-docs-blue.svg)](https://mpm.cb-geo.com/)
[![CircleCI](https://circleci.com/gh/cb-geo/mpm.svg?style=svg)](https://circleci.com/gh/cb-geo/mpm)
[![codecov](https://codecov.io/gh/cb-geo/mpm/branch/develop/graph/badge.svg)](https://codecov.io/gh/cb-geo/mpm)
[![](https://img.shields.io/github/issues-raw/cb-geo/mpm.svg)](https://github.com/cb-geo/mpm/issues)
[![Coverity](https://scan.coverity.com/projects/14389/badge.svg)](https://scan.coverity.com/projects/14389/badge.svg)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/cb-geo/mpm.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/cb-geo/mpm/context:cpp)
[![Project management](https://img.shields.io/badge/projects-view-ff69b4.svg)](https://github.com/cb-geo/mpm/projects/)

## Documentation

Please refer to [CB-Geo MPM Documentation](https://cb-geo.github.io/mpm-doc) for information on compiling, and running the code. The documentation also include the MPM theory.

## Install dependencies

* Docker image for CB-Geo mpm code [https://hub.docker.com/r/cbgeo/mpm](https://hub.docker.com/r/cbgeo/mpm)

* Instructions for running mpm docker container: [https://github.com/cb-geo/docker-mpm/blob/master/README.md](https://github.com/cb-geo/mpm-container/blob/master/README.md).

### Prerequisite packages
> The following prerequisite packages can be found in the docker image:

* [Boost](http://www.boost.org/)
* [Eigen](http://eigen.tuxfamily.org/)
* [Intel TBB](https://www.threadingbuildingblocks.org/)
* [HDF5](https://support.hdfgroup.org/HDF5/)

#### Optional
* [MPI](https://www.open-mpi.org/)
* [METIS](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/)
* [ParMETIS](http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/)
* [VTK](https://www.vtk.org/)

### Fedora installation

Please run the following command:

```shell
dnf install -y boost boost-devel clang cmake cppcheck eigen3-devel findutils gcc gcc-c++ \
                   git hdf5 hdf5-devel hdf5-openmpi hdf5-openmpi-devel kernel-devel lcov\
                   make openmpi openmpi-devel sqlite sqlite-devel tar tbb tbb-devel valgrind vim \
                   voro++ voro++-devel vtk vtk-devel wget
```

### Ubuntu installation

Please run the following commands to install dependencies:

```
sudo apt-get install -y gcc git libboost-all-dev libeigen3-dev libhdf5-serial-dev libopenmpi-dev \
                        libtbb-dev

```

To install other dependencies:
> CMake 3.15
```
sudo apt-get install software-properties-common
sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
sudo apt update
sudo apt upgrade
```

> OpenGL and X11:Xt
```
sudo apt-get install freeglut3-dev libxt-dev
```

> VTK
```
git clone git://vtk.org/VTK.git VTK
cd VTK && mkdir build && cd build/
cmake -DCMAKE_BUILD_TYPE:STRING=Release ..
make -j
sudo make install
```

### METIS/ParMETIS installation

```shell
# METIS and PARMETIS

wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz && \
    tar -xf metis-5.1.0.tar.gz && \
    cd metis-5.1.0/ && mkdir -p ~/workspace/metis && \
    make config shared=1 cc=mpicc cxx=mpic++ prefix=~/workspace/metis && \
    make install 

wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz && \
    tar -xf parmetis-4.0.3.tar.gz && \
    cd parmetis-4.0.3/ && mkdir -p ~/workspace/parmetis && \
    make config shared=1 cc=mpicc cxx=mpic++ prefix=~/workspace/parmetis && \
    make install
```



## Compile
> See https://mpm-doc.cb-geo.com/ for more detailed instructions. 

0. Run `mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ ..`.

1. Run `make clean && make -jN` (where N is the number of cores).

### Compile mpm or mpmtest

* To compile either `mpm` or `mpmtest` alone, run `make mpm -jN` or `make mpmtest -jN` (where N is the number of cores).

### Compile without tests [Editing CMake options]

To compile without tests run: `mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DMPM_BUILD_TESTING=Off  -DCMAKE_CXX_COMPILER=g++ -DMETIS_DIR=~/workspace/metis/ -DPARMETIS_DIR=~/workspace/parmetis/ ..`.

### Run tests

0. Run `./mpmtest -s` (for a verbose output) or `ctest -VV`.

### Run MPM
> See https://mpm-doc.cb-geo.com/ for more detailed instructions. 

The CB-Geo MPM code uses a `JSON` file for input configuration. To run the mpm code:

```
   ./mpm  [-p <tbb_parallel>] [-i <input_file>] -f <working_dir> [--]
          [--version] [-h]
```

For example:

```
./mpm -f /path/to/input-dir/ -i mpm-usf-3d.json -p 8
```

Where:

```

   -p <tbb_parallel>,  --tbb_parallel <tbb_parallel>
     Number of parallel TBB threads

   -i <input_file>,  --input_file <input_file>
     Input JSON file [mpm.json]

   -f <working_dir>,  --working_dir <working_dir>
     (required)  Current working folder

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

## Compile with MPI (Running on a cluster)

The CB-Geo mpm code can be compiled with `MPI` to distribute the workload across compute nodes in a cluster.

Additional steps to load `OpenMPI` on Fedora:

```
source /etc/profile.d/modules.sh
export MODULEPATH=$MODULEPATH:/usr/share/modulefiles
module load mpi/openmpi-x86_64
```

Compile with OpenMPI:

```
mkdir build && cd build 
export CXX_COMPILER=mpicxx
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ -DMETIS_DIR=~/workspace/metis/ -DPARMETIS_DIR=~/workspace/parmetis/ ..
make -jN
```

### Running the code with MPI

To run the CB-Geo mpm code on a cluster with MPI:

```
mpirun -N <#-MPI-tasks> ./mpm -f /path/to/input-dir/ -i mpm.json
```

For example to run the code on 4 compute nodes (MPI tasks):

```
mpirun -N 4 ./mpm -f ~/benchmarks/3d/uniaxial-stress -i mpm.json
```


## Citation

Kumar, K., Salmond, J., Kularathna, S., Wilkes, C., Tjung, E., Biscontin, G., & Soga, K. (2019). Scalable and modular material point method for large scale simulations. 2nd International Conference on the Material Point Method. Cambridge, UK.
