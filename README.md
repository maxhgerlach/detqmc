           /$$             /$$
          | $$            | $$
      /$$$$$$$  /$$$$$$  /$$$$$$    /$$$$$$  /$$$$$$/$$$$   /$$$$$$$
     /$$__  $$ /$$__  $$|_  $$_/   /$$__  $$| $$_  $$_  $$ /$$_____/
    | $$  | $$| $$$$$$$$  | $$    | $$  \ $$| $$ \ $$ \ $$| $$
    | $$  | $$| $$_____/  | $$ /$$| $$  | $$| $$ | $$ | $$| $$
    |  $$$$$$$|  $$$$$$$  |  $$$$/|  $$$$$$$| $$ | $$ | $$|  $$$$$$$
     \_______/ \_______/   \___/   \____  $$|__/ |__/ |__/ \_______/
                                        | $$
                                        | $$       generic C++ replica exchange
                                        |__/       determinantal quantum Monte Carlo code


generic nature of the code, allows to plug in a replica class
implementing any suitable model and to use that in numerically
stabilized finite temperature simulations with a replica exchange
mechanism

specifically implemented metallic O(N) spin-density wave (SDW) model,
N = 1,2,3; also includes some starter code for a Hubbard model replica

high - performance

list of what is implemented

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
    - [Single replica simulation](#single-replica-simulation)
    - [Replica exchange simulation](#replica-exchange-simulation)
- [Structure of DQMC code](#structure-of-dqmc-code)
- [Executables](#executables)
- [Overview of source code files](#overview-of-source-code-files)

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

## Single replica simulation ##

An example single replica simulation run is prepared in
`example/sdw-o2-single-replica/simulation.conf`.  From the project
root do the following to execute it:

``` shell
$ cd example/sdw-o2-single-replica
$ ../../Release/detqmcsdwopdim
```

This will set up the simulation according to the parameter values
specified in `simulation.conf`.  You can get an annotated list of all
possible command line or config file options with `detqmcsdwopdim
--help`.

Once the simulation has started, some information is printed to the
screen.  After it has finished, you will find some new files in the
same directory: `info.dat` summarizes data about the simulation.
`results.values` contains expectation values with error bars for some
measured observables.  The time series of these observable
measurements are contained in files `*.series`.  Then there is a file
`configs-phi.binarystream` containing the raw system configurations
sampled in the simulation to be evaluated subsequently.  Finally,
`simulation.state` contains all the checkpointing information to
resume a simulation that has been interrupted previously.  By running
`detqmcsdwopdim --sweeps 200` it is also possible to continue a
previously finished simulation until a higher target sweep count is
reached.

## Replica exchange simulation ##

A slightly more involved example for a replica exchange simulation
with eight different values of the tuning parameter r is prepared in
`example/sdw-o2-replica-exchange/simulation.conf`.  Typically, you
would run it like this:

``` shell
$ cd example/sdw-o2-replica-exchange
$ mpirun -n 8 ../../Release/detqmcptsdwopdim
```

Depending on your environment, you may need to adapt the invocation of
MPI, with the SLURM scheduler you would use `srun` instead of `mpirun`
for instance.

This will produce a common `info.dat` for all replicas, three files
`exchange-*.values` with statistics on the replica exchange process,
and a file `simulation.*.state` for each replica such that the
entire simulation can be resumed again.  Also, a number of
subdirectories has appeared, one for each value of the tuning
parameter r.  In these subdirectories the measurements corresponding
to that value of r are collected.

See `scripts/detptsubmit.py` for some inspiration how to set up a
number of large scale simulation jobs more easily, taking a file like
`example/simulation.job` as input.

# Structure of DQMC code #

Non-model specific aspects of the simulation are handled by class
`DetQMC` ([`detqmc.h`](src/detqmc.h)) for single replica simulations
and class `DetQMCPT` ([`detqmcpt.h`](src/detqmcpt.h)) for replica
exchange simulations.  They depend on a template parameter `Model`,
which should be filled with a class that implements the interface of
the abstract class `DetModel` ([`detmodel.h`](src/detmodel.h)).

Generally, you would want to use the class `DetModelGC`
([`detmodel.h`](src/detmodel.h)), derived from `DetModel`, and further
specialize that in a derived class.  It provides skeleton routines for
the numerically stable computation of the Green's function in a DQMC
sweep.  These skeleton routines are brought to life by the derived
class implementing a mdoel, which provides the routines that depend on
its fermion matrix structure.

The best example for such a class is `DetSDW`
([`detsdwopdim.h`](src/detsdwopdim.h)) for the metallic O(N)
spin-density wave (SDW) model with N = `OPDIM` = 1, 2, 3.  It also
implements the interface that's needed for the replica exchange scheme
managed by `DetQMCPT`.

A less complete and less optimized implementation of the half-filled
repulsive Hubbard model is also available in `DetHubbard`
([`dethubbard.h`](src/dethubbard.h)).

Observable values measured in a model are shared via instances of
`Observable` ([`observable.h`](src/observable.h)) such that the
observable handler classes
([`observablehandler.h`](src/observablehandler.h) and
[`mpiobservablehandlerpt.h`](src/mpiobservablehandlerpt.h)), managed
by `DetQMC[PT]`, can process them independently.

# Executables #

Most of these support an option `--help` to list supported command
line and configuration file options.

  * **Single replica DQMC simulations**
      * SDW model: `detqmcsdwopdim` supports order parameter
      dimensions N = 1, 2, and 3.  `detqmcsdwo1`, `detqmcsdwo2`, and
      `detqmcsdwo3` only support one of those, for the benefit of
      faster compilation.
      * Hubbard model: `detqmchubbard`
  * **Replica exchange (parallel tempering) DQMC simulations**
      * SDW model: `detqmcptsdwopdim`, `detqmcptsdwo1`,
        `detqmcptsdwo2`, `detqmcptsdwo3`
  * **Data evaluation**: `deteval`, `detevalbc`, `sdwcorr`,
    `sdweqtimesusc`, `tauintsimple`
  * **Multiple historgram reweighting**: `mrpt`, `mrptbc`,
    `mrpt_findIntersect`, `mrpt_binderRatioIntersect`,
    `mrptbc-binderratio-intersect`
  * **Utilities for data management**: `binarystreamtotext`,
    `binarystreamtonormmeanseries`,
    `binarystreamtonormmeanseriesrepeated`, `extractfrombinarystream`,
    `jointimeseries`


# Overview of source code files #
C++ source files are put in `src/`.  A number of Python scripts for
job management and data evaluation are collected in `scripts/`.

  * **Replica classes**, implementing one instance of a model in a Monte Carlo simulation
      * Generic base class for a model, implementing the skeleton of a
        numerically stabilized Monte Carlo sweep:
        [`detmodel.h`](src/detmodel.h),
        [`detmodel.cpp`](src/detmodel.cpp)
      * Generic model parameter handling:
        [`detmodelparams.h`](src/detmodelparams.h),
        [`detmodelloggingparams.h`](src/detmodelloggingparams.h),
        [`detmodelloggingparams.cpp`](src/detmodelloggingparams.cpp)
      * Helper struct: [`udv.h`](src/udv.h)
      * Metallic SDW model with O(1), O(2), or O(3) order parameter
        dimension: [`detsdwopdim.h`](src/detsdwopdim.h),
        [`detsdwopdim.cpp`](src/detsdwopdim.cpp),
        [`detsdwo1.cpp`](src/detsdwo1.cpp),
        [`detsdwo2.cpp`](src/detsdwo2.cpp),
        [`detsdwo3.cpp`](src/detsdwo3.cpp),
        [`detsdwparams.h`](src/detsdwparams.h),
        [`detsdwparams.cpp`](src/detsdwparams.cpp)
      * Serializing and passing around SDW model system
        configurations:
        [`detsdwsystemconfig.h`](src/detsdwsystemconfig.h),
        [`detsdwsystemconfig.cpp`](src/detsdwsystemconfig.cpp),
        [`detsdwsystemconfigfilehandle.h`](src/detsdwsystemconfigfilehandle.h)
      * Hubbard model: [`dethubbard.h`](src/dethubbard.h),
        [`dethubbard.cpp`](src/dethubbard.cpp),
        [`dethubbardparams.h`](src/dethubbardparams.h),
        [`dethubbardparams.cpp`](src/dethubbardparams.cpp)
  * **Handling of DQMC simulations**
      * Handling of a single-replica DQMC simulation (thermalization,
        production sweeps, measurements, saving of state) in class
        `DetQMC`: [`detqmc.h`](src/detqmc.h),
        [`detqmcparams.h`](src/detqmcparams.h),
        [`detqmcparams.cpp`](src/detqmcparams.cpp)
      * Handling of parallelized replica-exchange DQMC simulations in
        class `DetQMCPT` (`DetQMC` on a parallelized scale, with
        additional replica exchange moves):
        [`detqmcpt.h`](src/detqmcpt.h),
        [`detqmcptparams.h`](src/detqmcptparams.h),
        [`detqmcptparams.cpp`](src/detqmcptparams.cpp)
      * Observable measurements
          * Class shared between replica and observable handler
            classes: [`observable.h`](src/observable.h)
          * Classes to handle observable measurements in single replica
            simulations:
            [`observablehandler.h`](src/observablehandler.h),
            [`observablehandler.cpp`](src/observablehandler.cpp)
          * Classes to handle observable measurements in replica exchange
            simulations:
            [`mpiobservablehandlerpt.h`](src/mpiobservablehandlerpt.h),
            [`mpiobservablehandlerpt.cpp`](src/mpiobservablehandlerpt.cpp)
  * **Simulation main programs**
      * Single-replica simulations:
        [`maindetqmcsdwopdim.cpp`](src/maindetqmcsdwopdim.cpp),
        [`maindetqmcsdwo1.cpp`](src/maindetqmcsdwo1.cpp),
        [`maindetqmcsdwo2.cpp`](src/maindetqmcsdwo2.cpp),
        [`maindetqmcsdwo3.cpp`](src/maindetqmcsdwo3.cpp),
        [`maindetqmchubbard.cpp`](src/maindetqmchubbard.cpp)
      * Replica-exchange simulations:
        [`mpimaindetqmcptsdwopdim.cpp`](src/mpimaindetqmcptsdwopdim.cpp),
        [`mpimaindetqmcptsdwo1.cpp`](src/mpimaindetqmcptsdwo1.cpp),
        [`mpimaindetqmcptsdwo2.cpp`](src/mpimaindetqmcptsdwo2.cpp),
        [`mpimaindetqmcptsdwo3.cpp`](src/mpimaindetqmcptsdwo3.cpp)
  * `mrpt`: **Multiple histogram reweighting** for parallel tempering / replica exchange simulations
      * Core `mrpt` routines: [`mrpt.h`](src/mrpt.h),
        [`mrpt.cpp`](src/mrpt.cpp)
      * Core `mrpt` routines with jackknife error handling:
        [`mrpt-jk.h`](src/mrpt-j.h), [`mrpt-jk.cpp`](src/mrpt-jk.cpp)
      * Intersection points of observable measurement curves:
        [`mrpt-find-intersect.h`](src/mrpt-find-intersect.h),
        [`mrpt-find-intersect.cpp`](src/mrpt-find-intersect.cpp),
        [`mrpt-binderratio-intersect.h`](src/mrpt-binderratio-intersect.h)
        [`mrpt-binderratio-intersect.cpp`](src/mrpt-binderratio-intersect.cpp),
      * High-level interface with SWIG Python bindings:
        [`mrpt-highlevel.h`](src/mrpt-highlevel.h),
        [`mrpt-highlevel.cpp`](src/mrpt-highlevel.cpp),
        [`mrpt-highlevel.i`](src/mrpt-highlevel.i)
      * Main `mrpt` programs: [`main-mrpt.cpp`](src/main-mrpt.cpp),
        [`main-mrptbc.cpp`](src/main-mrptbc.cpp),
        [`main-mrpt-find-intersect.cpp`](src/main-mrpt-find-intersect.cpp),
        [`main-mrpt-binderratio-intersect.cpp`](src/main-mrpt-binderratio-intersect.cpp),
        [`main-mrptbc-binderratio-intersect.cpp`](src/main-mrptbc-binderratio-intersect.cpp)
      * Helper structs:
        [`reweightingresult.h`](src/reweightingresult.h),
        [`reweightedmomentsjk.h`](src/reweightedmomentsjk.h)
  * **Utility classes and functions**
      * Boost serialization support for additional classes:
      [`boost_serialize_armadillo.h`](src/boost_serialize_armadillo.h),
      [`boost_serialize_array.h`](src/boost_serialize_array.h),
      [`boost_serialize_uniqueptr.h`](src/boost_serialize_uniqueptr.h),
      [`boost_serialize_vector_uniqueptr.h`](src/boost_serialize_vector_uniqueptr.h)
      * `checkarray` derived from `std::array`, providing a bound
        checked `operator[]` in debug builds: [`checkarray.h`](src/checkarray.h)
      * Data IO to/from text files:
        [`datamapwriter.h`](src/datamapwriter.h),
        [`dataseriesloader.h`](src/dataseriesloader.h),
        [`dataserieswriter.h`](src/dataserieswriter.h),
        [`dataserieswritersucc.h`](src/dataserieswritersucc.h)
      * Exception handling: [`exceptions.h`](src/exceptions.h)
      * git revision and build information:
        [`git-revision.h`](src/git-revision.h),
        [`git-revision.c.in`](src/git-revision.c.in)
      * Histograms and related functions: [`histograms.h`](src/histograms.h)
      * Calcuations with internally logarithmic representation:
        [`logval.h`](src/logval.h)
      * Metadata for stored data: [`metadata.h`](src/metadata.h),
        [`metadata.cpp`](src/metadata.cpp)
      * Neighbor table for hypercubic lattices:
        [`neighbortable.h`](src/neighbortable.h)
      * Normal distributed random numbers:
        [`normaldistribution.h`](src/normaldistribution.h)
      * Routines for minimization and determination of roots:
        [`numerics.h`](src/numerics.h)
      * Run a linked Python interpreter to plot internal data:
        [`pytools.h`](src/pytools.h), [`pytools.cpp`](src/pytools.cpp)
      * Wrapper around random number generator:
        [`rngwrapper.h`](src/rngwrapper.h),
        [`rngwrapper.cpp`](src/rngwrapper.cpp)
      * Running averages: [`RunningAverage.h`](src/RunningAverage.h)
      * Statistics on time series: [`statistics.h`](src/statistics.h),
        [`statistics.cpp`](src/statistics.cpp)
      * Data structures for symmetric matrices:
        [`symmat.h`](src/symmat.h)
      * Optional execution timing of repeatedly executed code regions:
        [`timing.h`](src/timing.h), [`timing.cpp`](src/timing.cpp)
      * Various helpers: [`tools.h`](src/tools.h),
        [`tools.cpp`](src/tools.cpp)
      * Macros useful in `gdb` debugging sessions:
        [`toolsdebug.h`](src/toolsdebug.h)
      * Python-like `zip()` function for arbitrary containers: [`zip.h`](src/zip.h)
  * **Tools for data evaluation**
      * Expectation values from indivual time series:
        [`maindeteval.cpp`](src/maindeteval.cpp),
        [`maindetevalbc.cpp`](src/maindetevalbc.cpp)
      * Bosonic observables from stored system configurations:
        [`mainsdwcorr.cpp`](src/mainsdwcorr.cpp),
        [`mainsdweqtimesusc.cpp`](src/mainsdweqtimesusc.cpp)
      * Time series management:
        [`mainjointimeseries.cpp`](src/mainjointimeseries.cpp)
      * Autocorrelation times: [`maintauintsimple.cpp`](src/maintauintsimple.cpp)
      * Handling of stored streams of system configurations:
        [`mainbinarystreamtonormmeanseries.cpp`](src/mainbinarystreamtonormmeanseries.cpp),
        [`mainbinarystreamtonormmeanseriesrepeated.cpp`](src/mainbinarystreamtonormmeanseriesrepeated.cpp),
        [`mainbinarystreamtotext.cpp`](src/mainbinarystreamtotext.cpp),
        [`mainextractfrombinarystream.cpp`](src/mainextractfrombinarystream.cpp)
