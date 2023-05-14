cmake_minimum_required(VERSION 3.8.2)

project(QC VERSION 0.0.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
include(FetchContent)

include(CTest)

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
#Libint2
find_package(Libint2 REQUIRED)

add_subdirectory(src)
add_subdirectory(qc)
#add_subdirectory(tests)

# optionally add doxygen target to generate documentation
option(BUILD_DOCS "Enable building of documentation (requires Doxygen)" OFF)

if(BUILD_DOCS)
    find_package(Doxygen REQUIRED)
    set(DOXYGEN_EXCLUDE_PATTERNS "${CMAKE_SOURCE_DIR}/ext/*")
    doxygen_add_docs(doxygen
            ${CMAKE_SOURCE_DIR}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            )
endif()
