cmake_minimum_required(VERSION 3.10)

project(num_methods CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES main.cpp)

add_executable(num_methods ${SOURCES})
target_include_directories(num_methods PRIVATE include)