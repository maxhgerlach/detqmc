# detqmc #

## Building ##

### Prerequisites ###

This code has only been tested on Linux.  You will need the following
to get started:

  * A reasonably modern C++ compiler and standard library.  These are known to work well:
      * g++ 4.9.4
      * Clang 3.4
      * Intel 15.0 in conjunction with the libraries coming with a recent g++
      
    Earlier versions might do the job (Intel 13 should be fine; g++
    4.6 may compile with some tinkering to work around C++11
    shortcomings, while some intermediate versions < 4.9.4 fail due to
    a bug).
  * [CMake](https://cmake.org/) version 2.8.5+
  * Implementations of [BLAS](http://www.netlib.org/blas/) and [LAPACK](http://www.netlib.org/lapack/)
      * [OpenBlas](http://www.openblas.net/) and LAPACK packages work well
      * Intel [MKL](https://software.intel.com/en-us/intel-mkl) is a fast commercial alternative (commonly installed on compute clusters)
  * [FFTW](http://www.fftw.org/) 3.x
      * you can link to Intel MKL instead
  * For replica exchange simulations: an implementation of [MPI](http://www.mcs.anl.gov/research/projects/mpi/)
    * [Open MPI](https://www.open-mpi.org/) works well
    * Intel MPI is fine too

All of these should be available via your distribution's package manager or in the form of modules on your HPC cluster.
