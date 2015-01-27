# turn off warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

add_library(boost_system error_code.cpp local_free_on_destruction.hpp)