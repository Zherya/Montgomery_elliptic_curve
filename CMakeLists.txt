cmake_minimum_required(VERSION 3.8)
project(Montgomery_elliptic_curve)

include_directories(/usr/local/lib /usr/local/include)

set(CMAKE_CXX_STANDARD 17)

find_library(GMP_LIBRARY gmp)

set(SOURCE_FILES main.cpp Montgomery_curve.hpp Montgomery_curve.cpp)
add_executable(Montgomery_elliptic_curve ${SOURCE_FILES})

target_link_libraries(Montgomery_elliptic_curve "${GMP_LIBRARY}")