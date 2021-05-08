cmake_minimum_required(VERSION 3.16)
project(binary_row)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(src)

add_executable(binary_row_test example.cpp)
#target_link_libraries(binary_row_test binary_row_library)

add_subdirectory(tests)