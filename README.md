# 2D/3D Material Point Method (mpm)
> Cambridge Berkeley - Geomechanics

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/cb-geo/mpm/develop/license.md)
[![Developer docs](https://img.shields.io/badge/developer-docs-blue.svg)](http://cb-geo.github.io/mpm)
[![User docs](https://img.shields.io/badge/user-docs-blue.svg)](https://mpm.cb-geo.com/)
[![CircleCI](https://circleci.com/gh/cb-geo/mpm.svg?style=svg)](https://circleci.com/gh/cb-geo/mpm)
[![codecov](https://codecov.io/gh/cb-geo/mpm/branch/develop/graph/badge.svg)](https://codecov.io/gh/cb-geo/mpm)
[![](https://img.shields.io/github/issues-raw/cb-geo/mpm.svg)](https://github.com/cb-geo/mpm/issues)
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
* [VTK](https://www.vtk.org/)

### Fedora installation

Please run the following command:

```shell
dnf install -y boost boost-devel clang cmake cppcheck eigen3-devel findutils gcc gcc-c++ \
                   git hdf5 hdf5-devel kernel-devel lcov\
                   make openmpi openmpi-devel sqlite sqlite-devel tar tbb tbb-devel valgrind vim \
                   voro++ voro++-devel vtk vtk-devel wget
```

## Compile
> See https://mpm-doc.cb-geo.com/ for more detailed instructions. 

0. Run `mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release /path/to/CMakeLists.txt`.

1. Run `make clean && make -jN` (where N is the number of cores).

### Compile mpm or mpmtest

* To compile either `mpm` or `mpmtest` alone, run `make mpm -jN` or `make mpmtest -jN` (where N is the number of cores).

### Compile without tests [Editing CMake options]

To compile without tests run: `mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DMPM_BUILD_TESTING=Off /path/to/CMakeLists.txt`.


### Run tests

0. Run `./mpmtest -s` (for a verbose output) or `ctest -VV`.

### Run MPM
> See https://mpm-doc.cb-geo.com/ for more detailed instructions. 

The CB-Geo MPM code uses a `JSON` file for input configuration. To run the mpm code:

```
./mpm  -a <analysis> [-i <input_file>] -f <working_dir> [--] [--version]
       [-h]
```

Where:

```
   -a <analysis>,  --analysis <analysis>
     (required)  MPM analysis 
     (MPMExplicitUSF3D / MPMExplicitUSF2D / MPMExplicitUSL3D / MPMExplicitUSL2D)

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

## References
* [Aspect](https://github.com/geodynamics/aspect)
* [Dealii](https://github.com/dealii/dealii)
