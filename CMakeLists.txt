cmake_minimum_required(VERSION 3.17)
project(Task2)

set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CUDA_STANDARD 11)

include_directories(.)

add_executable(Task2
        src/functions.h
        src/Grid3D.cpp
        src/Grid3D.h
        src/utils.h
        src/main.cpp
        src/MathSolver.cpp
        src/MathSolver.h
        src/CpuSolver.cpp
        src/CpuSolver.h
        src/Block.cpp src/Block.h
        src/MPIProxy.cpp src/MPIProxy.h
        src/CudaSolver.cuh src/CudaSolver.cu)

#set_target_properties(
#        Task2
#        PROPERTIES
#        CUDA_SEPARABLE_COMPILATION ON)