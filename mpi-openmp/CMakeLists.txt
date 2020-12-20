cmake_minimum_required(VERSION 3.17)
project(Task2)

set(CMAKE_CXX_STANDARD 98)

include_directories(..)

add_executable(Task2
        src/functions.h
        src/Grid3D.cpp
        src/Grid3D.h
        src/utils.h
        src/main.cpp
        src/MathSolver.cpp
        src/MathSolver.h src/Block.cpp src/Block.h src/MPIProxy.cpp src/MPIProxy.h)
