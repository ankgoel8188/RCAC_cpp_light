cmake_minimum_required(VERSION 3.0)

project(RCAC)

message(STATUS "xxxxxxxxxxxxxxxxx  ${PROJECT_SOURCE_DIR}")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-Wall")
include_directories(src)

add_subdirectory(src)
