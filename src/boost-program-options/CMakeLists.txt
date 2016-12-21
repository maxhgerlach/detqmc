# turn off warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

add_library(boost_program_options cmdline.cpp config_file.cpp convert.cpp options_description.cpp parsers.cpp positional_options.cpp split.cpp utf8_codecvt_facet.cpp value_semantic.cpp variables_map.cpp winmain.cpp)