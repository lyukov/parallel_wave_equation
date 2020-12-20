cmake_minimum_required(VERSION 3.17)
project(Task2 CUDA)

set(CMAKE_CUDA_STANDARD 11)

add_executable(Task2 src/main.cu)

set_target_properties(
        Task2
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)