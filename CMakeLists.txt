cmake_minimum_required(VERSION 3.23)
project(P2)

set(CMAKE_CXX_STANDARD 17)

#Prerequisites
find_package(Eigen3 REQUIRED)
find_package(Libint2 REQUIRED)

add_executable(P2 main.cpp hf.cpp)

target_link_libraries(P2 Eigen3::Eigen)
target_link_libraries(P2 Libint2::cxx)