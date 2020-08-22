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
[![Project management](https://img.shields.io/badge/projects-view-ff69b4.svg)](https://github.com/orgs/cb-geo/projects/1)
[![Discourse forum](https://img.shields.io/badge/forum-mpm-blueviolet.svg)](https://cb-geo.discourse.group/c/mpm/)

## Documentation

Please refer to [CB-Geo MPM Documentation](https://mpm.cb-geo.com/) for information on compiling, and running the code. The documentation also include the MPM theory.

If you have any issues running or compiling the MPM code please open a issue on the [CB-Geo Discourse forum](https://cb-geo.discourse.group/c/mpm/).

## Running code on Docker

* Docker image for CB-Geo mpm code [https://hub.docker.com/r/cbgeo/mpm](https://hub.docker.com/r/cbgeo/mpm)

* Instructions for running mpm docker container: [https://github.com/cb-geo/docker-mpm/blob/master/README.md](https://github.com/cb-geo/mpm-container/blob/master/README.md).

## Running code locally

### Prerequisite packages
> The following prerequisite packages can be found in the docker image:

* [Boost](http://www.boost.org/)
* [Eigen](http://eigen.tuxfamily.org/)
* [HDF5](https://support.hdfgroup.org/HDF5/)

#### Optional
* [MKL](https://software.intel.com/en-us/mkl)
* [MPI](https://www.open-mpi.org/)
* [OpenMP 5.0](https://www.openmp.org/specifications/)
* [KaHIP](https://github.com/schulzchristian/KaHIP)
* [Partio](https://github.com/wdas/partio)
* [VTK](https://www.vtk.org/)

### Fedora installation (recommended)

Please run the following command:

```shell
dnf install -y boost boost-devel clang clang-analyzer clang-tools-extra cmake cppcheck dnf-plugins-core \
                   eigen3-devel findutils freeglut freeglut-devel gcc gcc-c++ git hdf5 hdf5-devel \
                   kernel-devel lcov libnsl make ninja-build openmpi openmpi-devel tar \
                   valgrind vim vtk vtk-devel wget
```

### Ubuntu installation

Please run the following commands to install dependencies:

```
sudo apt update
sudo apt upgrade
sudo apt install -y gcc git libboost-all-dev libeigen3-dev libhdf5-serial-dev libopenmpi-dev libomp-dev
```

If you are running Ubuntu 18.04 or below, you may want to update the GCC version to 9 to have OpenMP 5 specifications
support.

```
sudo apt install software-properties-common
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt install gcc-9 g++-9
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 90 --slave /usr/bin/g++ g++ /usr/bin/g++-9 --slave /usr/bin/gcov gcov /usr/bin/gcov-9

```

To install other dependencies:
> CMake 3.15
```
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
git clone https://gitlab.kitware.com/vtk/vtk.git VTK
cd VTK && mkdir build && cd build/
cmake -DCMAKE_BUILD_TYPE:STRING=Release ..
make -j
sudo make install
```

### Partio for Houdini SFX Visualization

```shell
mkdir -p ~/workspace && cd ~/workspace/ && git clone https://github.com/wdas/partio.git && \
    cd partio && cmake . && make
```

Houdini supported (*.bgeo) files will be generated. These can be rendered using the non-commercial [Houdini Apprentice](https://www.sidefx.com/download/).

### KaHIP installation for domain decomposition

```shell
cd ~/workspace/ && git clone https://github.com/schulzchristian/KaHIP && \
   cd KaHIP && sh ./compile_withcmake.sh
```

## Compile
> See [CB-Geo MPM Documentation](https://mpm.cb-geo.com/) for more detailed instructions.

0. Run `mkdir build && cd build && cmake -DCMAKE_CXX_COMPILER=g++ ..`.

1. Run `make clean && make -jN` (where N is the number of cores).

> To compile without KaHIP partitioning use `cmake -DNO_KAHIP=True ..`

### Compile mpm or mpmtest

* To compile either `mpm` or `mpmtest` alone, run `make mpm -jN` or `make mpmtest -jN` (where N is the number of cores).

### Compile without tests [Editing CMake options]

To compile without tests run: `mkdir build && cd build && cmake -DMPM_BUILD_TESTING=Off  -DCMAKE_CXX_COMPILER=g++ ..`.

## Compile with MPI (Running on a cluster)

The CB-Geo mpm code can be compiled with `MPI` to distribute the workload across compute nodes in a cluster.

Additional steps to load `OpenMPI` on Fedora:

```
source /etc/profile.d/modules.sh
export MODULEPATH=$MODULEPATH:/usr/share/modulefiles
module load mpi/openmpi-x86_64
```

Compile with OpenMPI (with halo exchange):

```
mkdir build && cd build
export CXX_COMPILER=mpicxx
cmake -DCMAKE_BUILD_TYPE=Release -DKAHIP_ROOT=~/workspace/KaHIP/ -DHALO_EXCHANGE=On ..
make -jN
```

To enable halo exchange set `-DHALO_EXCHANGE=On` in `CMake`. Halo exchange is a better MPI communication protocol, however, use this only for larger number of MPI tasks (> 4).

### Compile with Ninja build system [Alternative to Make]

0. Run `mkdir build && cd build && cmake -GNinja -DCMAKE_CXX_COMPILER=g++ ..`.

1. Run `ninja`

### Compile with Partio viz support

Please include `-DPARTIO_ROOT=/path/to/partio/` in the cmake command. A typical cmake command would look like `cmake -DCMAKE_BUILD_TYPE=Release -DPARTIO_ROOT=~/workspace/partio/ ..`

## Run tests

0. Run `./mpmtest -s` (for a verbose output) or `ctest -VV`.

## Run MPM
> See [CB-Geo MPM Documentation](https://mpm.cb-geo.com/) for more detailed instructions.

The CB-Geo MPM code uses a `JSON` file for input configuration. To run the mpm code:

```
   ./mpm  [-p <parallel>] [-i <input_file>] -f <working_dir> [--]
          [--version] [-h]
```

For example:

```
export OMP_SCHEDULE="static,4"
./mpm -f /path/to/input-dir/ -i mpm-usf-3d.json
```

Where:

```
   -p <parallel>,  --parallel <parallel>
     Number of parallel threads

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

### Running the code with MPI

To run the CB-Geo mpm code on a cluster with MPI:

```
mpirun -N <#-MPI-tasks> ./mpm -f /path/to/input-dir/ -i mpm.json
```

For example to run the code on 4 compute nodes (MPI tasks):

```
mpirun -N 4 ./mpm -f ~/benchmarks/3d/uniaxial-stress -i mpm.json
```

## Authors

Please refer to the [list of contributors to the CB-Geo MPM code](AUTHORS.md).

## Citation

If you publish results using our code, please acknowledge our work by quoting the following paper:

Kumar, K., Salmond, J., Kularathna, S., Wilkes, C., Tjung, E., Biscontin, G., & Soga, K. (2019). Scalable and modular material point method for large scale simulations. 2nd International Conference on the Material Point Method. Cambridge, UK. [https://arxiv.org/abs/1909.13380](https://arxiv.org/abs/1909.13380)
