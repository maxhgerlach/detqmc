detqmc title

overview

generic nature of the code, allows to plug in a replica class
implementing any suitable model and to use that in numerically
stabilized finite temperature simulations with a replica exchange
mechanism

specifically implemented metallic O(N) spin-density wave (SDW) model,
N = 1,2,3; also includes some starter code for a Hubbard model replica

refer to papers, thesis

include in overview: mention multiple histogram reweighting for data
analysis

MPL 2.0

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-generate-toc again -->
**Table of Contents**

- [Setup](#setup)
    - [Prerequisites](#prerequisites)
    - [Included libraries](#included-libraries)
    - [Compilation](#compilation)
- [Usage](#usage)
- [Overview over source code files](#overview-over-source-code-files)

<!-- markdown-toc end -->


# Setup #

## Prerequisites ##

This code has only been tested on Linux.  You will need the following
to get started:

  * A reasonably modern C++ compiler and standard library.  These are known to work well:
      * g++ 4.9.4
      * Clang 3.4
      * Intel 15.0 in conjunction with the libraries coming with a recent g++
      
    Earlier versions might do the job (Intel 13 should be fine; g++
    4.6 may compile with some tinkering to work around its C++11
    shortcomings, while some intermediate versions < 4.9.4 fail due to
    a bug).
  * [CMake](https://cmake.org/) version 2.8.5+
  * Implementations of [BLAS](http://www.netlib.org/blas/) and
    [LAPACK](http://www.netlib.org/lapack/)
      * [OpenBlas](http://www.openblas.net/) and LAPACK packages work well
      * Intel [MKL](https://software.intel.com/en-us/intel-mkl) is a
        fast commercial alternative (commonly installed on compute
        clusters)
  * [FFTW](http://www.fftw.org/) 3.x
      * you can link to Intel MKL instead
  * For replica exchange simulations: an implementation of
    [MPI](http://www.mcs.anl.gov/research/projects/mpi/)
    * [Open MPI](https://www.open-mpi.org/) works well
    * Intel MPI is fine too

All of these should be available via your distribution's package
manager or in the form of modules on your HPC cluster.

## Included libraries ##

For simplicity (portions of) these libraries are included with this package:

  * [Armadillo](http://arma.sourceforge.net/) 7.600.2 for high-level linear algebra
  * [dSFMT](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/#dSFMT) 2.1
    for pseudo-random number generation
  * [Boost](http://www.boost.org/) 1.51
  * [Dlib](http://dlib.net/), used in the mrpt reweighting code
  * [cnpy](https://github.com/rogersce/cnpy), used in sdwcorr

License information is included within the respective subdirectories.

## Compilation ##

To compile an optimized full build with default settings do:

``` shell
cd Release
./runcmake.sh
make -j
```

After having run CMake you can alternatively build only select
targets.  For example, run `make -j detqmcptsdwo2` to build all
requirements for the executable to run parallelized replica exchange
DQMC simulations of the O(2) SDW model.

Go to the directory `Debug` instead of `Release` for a build with
debug symbols, additional error checks, and no optimization.

If this does not work immediately, some tweaking may be necessary.

You can change the invocation  of `cmake` in `Release/runcmake.sh`.  For
instance, to use `clang` instead of `g++`, edit it to read

``` shell
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DUSE_GNU_COMPILER=OFF -DUSE_CLANG_COMPILER=ON  ../src
```
There are a couple more options defined in `src/CMakeList.txt`.

For more specific configuration you will have to directly edit the
`src/CMakeList.txt` file.  Its inital section contains logic for
system specific build settings.  It contains some example sections for
the "Cheops" and "Jureca" HPC clusters, and for some machines running
outdated system installations: the "l71" work station and the
computers in a domain "thp".  The parts of the CMake script that
enable these specific settings are currently commented out.  Use them
as example guidelines to adapt the setup to your own machines.

Depending on your BLAS and LAPACK installation you might also have to
uncomment some preprocessor directives in the file
`src/armadillo/armadillo_bits/config.hpp`.  Alternatively, you could
just delete the entire directory `src/armadillo` and use your own,
independent installation of [Armadillo](http://arma.sourceforge.net/).


# Usage #


# Overview over source code files #
  * Utility classes and functions
      * Boost serialization support for additional classes:
      [`boost_serialize_armadillo.h`](src/boost_serialize_armadillo.h),
      [`boost_serialize_array.h`](src/boost_serialize_array.h),
      [`boost_serialize_uniqueptr.h`](src/boost_serialize_uniqueptr.h),
      [`boost_serialize_vector_uniqueptr.h`](src/boost_serialize_vector_uniqueptr.h)
      * `checkarray` derived from `std::array`, providing a bound
        checked `operator[]` in debug builds: [`checkarray.h`](src/checkarray.h)
      * data IO to/from text files:
        [`datamapwriter.h`](src/datamapwriter.h),
        [`dataseriesloader.h`](src/dataseriesloader.h),
        [`dataserieswriter.h`](src/dataserieswriter.h),
        [`dataserieswritersucc.h`](src/dataserieswritersucc.h)
  * Replica classes, implementing one instance of a model in a Monte Carlo simulation
      * Generic base class for a model, implementing the skeleton of a
        numerically stabilized Monte Carlo sweep:
        [`detmodel.h`](/src/detmodel.h),
        [`detmodel.cpp`](/src/detmodel.cpp)
      * Generic model parameter handling:
        [`detmodelparams.h`](/src/detmodelparams.h),
        [`detmodelloggingparams.h`](/src/detmodelloggingparams.h),
        [`detmodelloggingparams.cpp`](/src/detmodelloggingparams.cpp)
      * Hubbard model: [`dethubbard.h`](src/dethubbard.h),
        [`dethubbard.cpp`](src/dethubbard.cpp),
        [`dethubbardparams.h`](src/dethubbardparams.h),
        [`dethubbardparams.cpp`](src/dethubbardparams.cpp)
