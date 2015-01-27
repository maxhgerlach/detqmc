# turn off warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

add_library(boost_filesystem codecvt_error_category.cpp operations.cpp path.cpp path_traits.cpp portability.cpp unique_path.cpp utf8_codecvt_facet.cpp windows_file_codecvt.cpp windows_file_codecvt.hpp)