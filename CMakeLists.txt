cmake_minimum_required(VERSION 3.22)
project(MEKF C)

find_package(Catch2 3 REQUIRED)

# set(CMAKE_C_STANDARD 11)
# set(CMAKE_EXE_LINKER_FLAGS "-static")
# set(CMAKE_C_COMPILER arm-linux-gnueabihf-gcc)

add_subdirectory(mesch12b)

add_executable(test test.c)