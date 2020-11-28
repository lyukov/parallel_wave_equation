cmake_minimum_required(VERSION 3.17)
project(Task2)

set(CMAKE_CXX_STANDARD 98)

include_directories(.)

add_executable(Task2
        functions.h
        Grid3D.cpp
        Grid3D.h
        log.h
        main.cpp
        Solver.cpp
        Solver.h Block.cpp Block.h)
