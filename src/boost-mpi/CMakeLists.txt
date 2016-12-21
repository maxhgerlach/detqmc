# turn off warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

add_library(boost_mpi
  broadcast.cpp communicator.cpp computation_tree.cpp content_oarchive.cpp environment.cpp exception.cpp graph_communicator.cpp group.cpp intercommunicator.cpp mpi_datatype_cache.cpp mpi_datatype_oarchive.cpp packed_iarchive.cpp packed_oarchive.cpp packed_skeleton_iarchive.cpp packed_skeleton_oarchive.cpp point_to_point.cpp request.cpp text_skeleton_oarchive.cpp timer.cpp)

# MPI compile / link flags
set_target_properties(boost_mpi PROPERTIES
  COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}"
  LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
