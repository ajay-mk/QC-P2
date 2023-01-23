cmake_minimum_required(VERSION 3.23)
project(P2)

set(CMAKE_CXX_STANDARD 17)
include(FetchContent)

#Prerequisites
# BTAS
FetchContent_Declare(
        BTAS
        GIT_REPOSITORY https://github.com/ValeevGroup/btas.git
)
FetchContent_MakeAvailable(BTAS)
FetchContent_GetProperties(BTAS
        SOURCE_DIR BTAS_SOURCE_DIR
        BINARY_DIR BTAS_BINARY_DIR
        )
# JSON
# Fetching release copy of nlohmann/json because of the large size of the original repo
FetchContent_Declare(json
        GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
        GIT_PROGRESS TRUE
        GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(json)

#find_package(Eigen3 REQUIRED)
find_package(Libint2 REQUIRED)


add_executable(P2 main.cpp hf.cpp mp2.cpp general.cpp general.h hf.h mp2.h integrals.cpp integrals.h cc.cpp cc.h)

#target_link_libraries(P2 Eigen3::Eigen)
target_link_libraries(P2 BTAS)
target_link_libraries(P2 nlohmann_json::nlohmann_json)
target_link_libraries(P2 Libint2::cxx)
