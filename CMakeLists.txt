cmake_minimum_required (VERSION 2.8.5)

# NOTE: For simplicity, this cmake script contains host-specific code
# and needs to be adapted for new build hosts.  In principle, it would
# be cleaner to adhere better to the CMake spirit and have these
# definitions outside this script (e.g. in a calling shell script).

# host discrimination
site_name(HOST_NAME)
# some specific settings for the Cheops cluster.
# this sets CHEOPS_POSITION to -1, if not found:
string(FIND "${HOST_NAME}" "cheops" CHEOPS_POSITION)
if ("${CHEOPS_POSITION}" LESS 0)
  set(CHEOPS_FOUND "FALSE")
else ("${CHEOPS_POSITION}" LESS 0)
  set(CHEOPS_FOUND "TRUE")
endif ("${CHEOPS_POSITION}" LESS 0)
# some specific settings for the l71 work station:
string(FIND "${HOST_NAME}" "l71" L71_POSITION)
if ("${L71_POSITION}" LESS 0)
  set(L71_FOUND "FALSE")
else ("${L71_POSITION}" LESS 0)
  set(L71_FOUND "TRUE")
endif ("${L71_POSITION}" LESS 0)
# some specific settings for the Thinkpad notebook:
string(FIND "${HOST_NAME}" "max-ThinkPad-X220" THINKPAD_POSITION)
if ("${THINKPAD_POSITION}" LESS 0)
  set(THINKPAD_FOUND "FALSE")
else ("${THINKPAD_POSITION}" LESS 0)
  set(THINKPAD_FOUND "TRUE")
endif ("${THINKPAD_POSITION}" LESS 0)

execute_process(COMMAND dnsdomainname OUTPUT_VARIABLE DOMAIN_NAME)
string(FIND "${DOMAIN_NAME}" "thp" THP_POSITION)
if ("${THP_POSITION}" LESS 0)
  set(THP_DOMAIN "FALSE")
else ("${THP_POSITION}" LESS 0)
  set(THP_DOMAIN "TRUE")
endif ("${THP_POSITION}" LESS 0)



# g++ 4.8 and 4.9 currently still have a bug preventing them to
# compile this code. So on the Thinkpad just use g++-4.7 manually
if (${THINKPAD_FOUND})
  set(CMAKE_C_COMPILER "gcc-4.7")
  set(CMAKE_CXX_COMPILER "g++-4.7")
endif (${THINKPAD_FOUND})



project (Det)


# standard CPPFLAGS etc


include_directories(
  "${PROJECT_SOURCE_DIR}"
  "$ENV{HOME}/include"
  "$ENV{HOME}/opt/include"
  )


# MPI
find_package(MPI)
# this will set (lang = C, CXX)
# MPI_<lang>_FOUND           TRUE if FindMPI found MPI flags for <lang>
# MPI_<lang>_COMPILER        MPI Compiler wrapper for <lang>
# MPI_<lang>_COMPILE_FLAGS   Compilation flags for MPI programs
# MPI_<lang>_INCLUDE_PATH    Include path(s) for MPI header
# MPI_<lang>_LINK_FLAGS      Linking flags for MPI programs
# MPI_<lang>_LIBRARIES       All libraries to link MPI programs against
if (NOT ${MPI_CXX_FOUND})
  message("No MPI CXX compiler found!")
else (NOT ${MPI_CXX_FOUND})
  include_directories(MPI_CXX_INCLUDE_PATH)
endif (NOT ${MPI_CXX_FOUND})
# -> remember to link ${MPI_CXX_LIBRARIES} for the respective
#    executables
# -> remember to set additional flags ${MPI_CXX_COMPILE_FLAGS},
#    ${MPI_CXX_LINK_FLAGS} for these targets



# Preprocessor flags: -D definitions
add_definitions(-DDSFMT -DDSFMT_MEXP=19937)
add_definitions(-DDUMA_NO_DUMA) # TODO SOMETIME reenable

# CXX flags:
# -march settings
if (NOT ${THP_DOMAIN})
  MARCH_FLAGS = "-march=native"
else (NOT ${THP_DOMAIN})
  # setting to use instead of --march=native -- to be able to run on Core2 computers in the thp domain
  # use a safer setting
  # [obtained via `cc -march=native -E -v - </dev/null 2>&1 | grep cc1` on l07]
  MARCH_FLAGS = "-march=core2 -mcx16 -msahf -mno-movbe -mno-aes -mno-pclmul -mno-popcnt -mno-abm -mno-lwp -mno-fma -mno-fma4 -mno-xop -mno-bmi -mno-tbm -mno-avx -mno-sse4.2 -msse4.1 --param l1-cache-size=32 --param l1-cache-line-size=64 --param l2-cache-size=6144 -mtune=core2 -fstack-protector"
endif (NOT ${THP_DOMAIN})
WARNING_FLAGS = "-Wall -Winit-self -Wextra -Woverloaded-virtual -Wconversion"
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${MARCH_FLAGS} ${WARNING_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MARCH_FLAGS} ${WARNING_FLAGS} -std=c++0x")

# build type specific flags:
# these are appended to ${CMAKE_C[XX]_FLAGS} for DEBUG / RELEASE / ETC builds

# TODO: support build configuration RelWithDebInfo
# TODO: support option to enable timing

# TODO: make sure only Intel compiler used on Cheops
if (${CHEOPS_FOUND})
  set(EXTRA_FLAGS_DEBUG "-g -O0 -fp-trap-common")
  # floating-point traps for division by zero, invalid operation, and overflow
else (${CHEOPS_FOUND})
  set(EXTRA_FLAGS_DEBUG "-g -O0")
endif (${CHEOPS_FOUND})
set(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   ${EXTRA_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${EXTRA_FLAGS_DEBUG}")
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  # TODO SOMETIME: option for duma
  #  link_directories("$ENV{HOME}/packages/duma_2_5_15")
  #  include_directories("$ENV{HOME}/packages/duma_2_5_15")
  #  set(EXTRA_LIBRARIES ${EXTRA_LIBRARIES} duma)
  add_definitions(-D_GLIBCXX_DEBUG -DMAX_DEBUG)
endif(CMAKE_BUILD_TYPE MATCHES "Debug")

# TODO: make sure only Intel compiler used on Cheops
if (${CHEOPS_FOUND})
  set(EXTRA_FLAGS_RELEASE "-O3 -ansi-alias -xSSE4.2 -no-prec-div")
else (${CHEOPS_FOUND})
  set(EXTRA_FLAGS_RELEASE "-O3 -fno-strict-aliasing -ffast-math")
endif (${CHEOPS_FOUND})
set(CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}   ${EXTRA_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_FLAGS_RELEASE}")
if(CMAKE_BUILD_TYPE MATCHES "Release")
  add_definitions(-DBOOST_DISABLE_ASSERTS -DARMA_NO_DEBUG -DNDEBUG -DHAVE_SSE2)
endif(CMAKE_BUILD_TYPE MATCHES "Release")



# External libraries -- Armadillo
set(EXTRA_LIBRARY_LOCATIONS "$ENV{HOME}/lib" "$ENV{HOME}/opt/lib")
if (${CHEOPS_FOUND})
  set(ARMADILLO_LIBRARIES
    "-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -liomp5 -lpthread")
  # do not link -lgfortran
  set(EXTRA_LIBRARY_LOCATIONS "${EXTRA_LIBRARY_LOCATIONS}" "/opt/rrzk/compiler/gnu/4.6.4/lib64")
else (${CHEOPS_FOUND})
  if (${L71_FOUND})
    set(ARMADILLO_LIBRARIES
      "-llapack -lblas -lgfortran")
  else (${L71_FOUND})
    set(ARMADILLO_LIBRARIES "-larmadillo")
  endif (${L71_FOUND})
  if (${THP_DOMAIN})
    # on thp machines: explicitly link to Lapack, not the Atlas implementation (which is not compatible to Openblas)
    set(EXTRA_LIBRARY_LOCATIONS "${EXTRA_LIBRARY_LOCATIONS}" "/usr/lib/lapack/")
  else (${THP_DOMAIN})
  endif (${THP_DOMAIN})
endif (${CHEOPS_FOUND})
link_directories(${EXTRA_LIBRARY_LOCATIONS})



# For machines that have this installed and easily accessible: include
# pytools.cpp and link to -lpython2.7
option(USE_PYTOOLS
  "include some python routines for interactivity (does not work on Cheops cluster)" ON)
if (USE_PYTOOLS)
  set(PYTOOLS_SRC "pytools.cpp")
  set(PYTHON_LIB "-lpython2.7")
else (USE_PYTOOLS)
  set(PYTOOLS_SRC "")
  set(PYTHON_LIB "")
endif (USE_PYTOOLS)


## foreign code in subdirectories -> have a separate CMakeLists.txt for those
# directly include needed boost source files to reduce linker woes
# (debug-builts don't work unless exactly the same compiler options
# are applied to the libraries linked in...)
add_subdirectory(boost-chrono)
add_subdirectory(boost-filesystem)
add_subdirectory(boost-mpi)
add_subdirectory(boost-program-options)
add_subdirectory(boost-serialization)
add_subdirectory(boost-system)
add_subdirectory(boost-timer)
# collect all boost libraries, but keep MPI separete:
set(BOOST_LIBS
  boost_program_options boost_filesystem boost_serialization
  boost_system boost_timer boost_chrono)
# we need -lrt for timer / chrono
set(BOOST_LIBS
  ${BOOST_LIBS} "-lrt")
add_subdirectory(dsfmt)


## Git revision info
# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Get the latest commit hash of the working branch
execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_REVISION_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)


configure_file(
  git-revision.c.in             # in source tree
  git-revision.c                # in binary tree - is this then found automatically?
)




## Executables

set(detqmchubbard_SRC
  dethubbard.cpp dethubbardparams.cpp
  detqmcparams.cpp maindetqmchubbard.cpp rngwrapper.cpp
  observablehandler.cpp metadata.cpp detmodel.cpp timing.cpp
  detmodelloggingparams.cpp)
set(detqmcsdwopdim_SRC
  maindetqmcsdwopdim.cpp rngwrapper.cpp observablehandler.cpp
  metadata.cpp detmodel.cpp detsdwopdim.cpp detsdwparams.cpp
  detmodelloggingparams.cpp detqmcparams.cpp timing.cpp)
set(detqmcsdwo1_SRC
  maindetqmcsdwo1.cpp rngwrapper.cpp
  observablehandler.cpp metadata.cpp detmodel.cpp detsdwo1.cpp
  detsdwparams.cpp detmodelloggingparams.cpp detqmcparams.cpp
  timing.cpp)
set(detqmcsdwo2_SRC
  maindetqmcsdwo2.cpp rngwrapper.cpp
  observablehandler.cpp metadata.cpp detmodel.cpp detsdwo2.cpp
  detsdwparams.cpp detmodelloggingparams.cpp detqmcparams.cpp
  timing.cpp)
set(detqmcsdwo3_SRC
  maindetqmcsdwo3.cpp rngwrapper.cpp
  observablehandler.cpp metadata.cpp detmodel.cpp detsdwo3.cpp
  detsdwparams.cpp detmodelloggingparams.cpp detqmcparams.cpp
  timing.cpp)
set(detqmcptsdwopdim_SRC
  mpimaindetqmcptsdwopdim.cpp
  rngwrapper.cpp mpiobservablehandlerpt.cpp metadata.cpp detmodel.cpp
  detsdwopdim.cpp detsdwparams.cpp detmodelloggingparams.cpp
  detqmcparams.cpp detqmcptparams.cpp timing.cpp)
set(detqmcptsdwo1_SRC
  mpimaindetqmcptsdwo1.cpp rngwrapper.cpp mpiobservablehandlerpt.cpp
  metadata.cpp detmodel.cpp detsdwo1.cpp detsdwparams.cpp
  detmodelloggingparams.cpp detqmcparams.cpp detqmcptparams.cpp
  timing.cpp)
set(detqmcptsdwo2_SRC
  mpimaindetqmcptsdwo2.cpp
  rngwrapper.cpp mpiobservablehandlerpt.cpp metadata.cpp detmodel.cpp
  detsdwo2.cpp detsdwparams.cpp detmodelloggingparams.cpp
  detqmcparams.cpp detqmcptparams.cpp timing.cpp)
set(detqmcptsdwo3_SRC
  mpimaindetqmcptsdwo3.cpp rngwrapper.cpp mpiobservablehandlerpt.cpp
  metadata.cpp detmodel.cpp detsdwo3.cpp detsdwparams.cpp
  detmodelloggingparams.cpp detqmcparams.cpp detqmcptparams.cpp
  timing.cpp)
set(binarystreamtotext_SRC
  mainbinarystreamtotext.cpp)
set(tauintsimple_SRC
  maintauintsimple.cpp metadata.cpp)
set(deteval_SRC
  maindeteval.cpp tools.cpp metadata.cpp)
set(detevalbc_SRC
  maindetevalbc.cpp tools.cpp metadata.cpp)
set(testspheresampling_SRC
  maintestspheresampling.cpp rngwrapper.cpp)

add_executable(detqmchubbard
  git-revision.c
  ${detqmchubbard_SRC}
  ${PYTOOLS_SRC}
  )
target_link_libraries(detqmchubbard
  dsfmt ${ARMADILLO_LIBRARIES} ${BOOST_LIBS} ${PYTHON_LIB})


add_executable(detqmcsdwopdim
  git-revision.c
  ${detqmcsdwopdim_SRC}
  ${PYTOOLS_SRC}
  )
target_link_libraries(detqmcsdwopdim
  dsfmt ${ARMADILLO_LIBRARIES} ${BOOST_LIBS} ${PYTHON_LIB})

add_executable(detqmcsdwo1
  git-revision.c
  ${detqmcsdwo1_SRC}
  ${PYTOOLS_SRC}
  )
target_link_libraries(detqmcsdwo1
  dsfmt ${ARMADILLO_LIBRARIES} ${BOOST_LIBS} ${PYTHON_LIB})

add_executable(detqmcsdwo2
  git-revision.c
  ${detqmcsdwo2_SRC}
  ${PYTOOLS_SRC}
  )
target_link_libraries(detqmcsdwo2
  dsfmt ${ARMADILLO_LIBRARIES} ${BOOST_LIBS} ${PYTHON_LIB})

add_executable(detqmcsdwo3
  git-revision.c
  ${detqmcsdwo3_SRC}
  ${PYTOOLS_SRC}
  )
target_link_libraries(detqmcsdwo3
  dsfmt ${ARMADILLO_LIBRARIES} ${BOOST_LIBS} ${PYTHON_LIB})


if ("${MPI_CXX_FOUND}")
  set(PT_SIM_LIBS
    ${MPI_CXX_LIBRARIES} dsfmt ${ARMADILLO_LIBRARIES}
    ${BOOST_LIBS} boost_mpi ${PYTHON_LIB})
  
  add_executable(detqmcptsdwopdim
    git-revision.c
    ${detqmcptsdwopdim_SRC}
    ${PYTOOLS_SRC}
    )
  target_link_libraries(detqmcptsdwopdim
    ${PT_SIM_LIBRARIES})

  add_executable(detqmcptsdwo1
    git-revision.c
    ${detqmcptsdwo1_SRC}
    ${PYTOOLS_SRC}
    )
  target_link_libraries(detqmcptsdwo1
    ${PT_SIM_LIBRARIES})

  add_executable(detqmcptsdwo2
    git-revision.c
    ${detqmcptsdwo2_SRC}
    ${PYTOOLS_SRC}
    )
  target_link_libraries(detqmcptsdwo2
    ${PT_SIM_LIBRARIES})

  add_executable(detqmcptsdwo3
    git-revision.c
    ${detqmcptsdwo3_SRC}
    ${PYTOOLS_SRC}
    )
  target_link_libraries(detqmcptsdwo3
    ${PT_SIM_LIBRARIES})

  # Compile / Link flags
  set(MPI_EXECUTABLES
    detqmcptsdwopdim detqmcptsdwo1 detqmcptsdwo2 detqmcptsdwo3)
  set_target_properties(${MPI_EXECUTABLES} PROPERTIES
    COMPILE_FLAGS ${MPI_CXX_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_CXX_LINK_FLAGS})
endif ("${MPI_CXX_FOUND}")

add_executable(binarystreamtotext
  ${binarystreamtotext_SRC})

add_executable(tauintsimple
  ${tauintsimple_SRC})

add_executable(deteval
  git-revision.c
  ${deteval_SRC}
  )
target_link_libraries(deteval
  boost_program_options)

add_executable(detevalbc
  git-revision.c
  ${detevalbc_SRC}
  )
target_link_libraries(detevalbc
  boost_program_options boost_filesystem boost_system)

add_executable(testspheresampling
  ${testspheresampling_SRC}
  )
taget_link_libraries(testspheresampling
  ${ARMADILLO_LIBRARIES} dsfmt)