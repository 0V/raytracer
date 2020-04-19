cmake_minimum_required (VERSION 3.1)

project(raytracer)

# For debug
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(BUILD_TARGET "raytracer")

set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++17")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_MINSIZEREL, "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE, "-O4 -DNDEBUG")

file(GLOB "${BUILD_TARGET}_SOURCES" ${PROJECT_SOURCE_DIR}/src/main.cpp)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/src)

include_directories(${BUILD_TARGET} "${INCLUDE_DIR}")

add_executable(${BUILD_TARGET} ${${BUILD_TARGET}_SOURCES})


