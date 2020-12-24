#include "CudaSolver.cuh"
#include "utils.h"
#include <cstdlib>
#include <cuda.h>
#include <cmath>

#define SAFE_CALL(CallInstruction) { \
    cudaError_t cuerr = CallInstruction; \
    if(cuerr != cudaSuccess) { \
        printf("CUDA error: %s at call \"" #CallInstruction "\"\n", cudaGetErrorString(cuerr)); \
            throw "error in CUDA API function, aborting..."; \
    } \
}

#define SAFE_KERNEL_CALL(KernelCallInstruction) { \
    KernelCallInstruction; \
    cudaError_t cuerr = cudaGetLastError(); \
    if(cuerr != cudaSuccess) { \
        printf("CUDA error in kernel launch: %s at kernel \"" #KernelCallInstruction "\"\n", cudaGetErrorString(cuerr)); \
            throw "error in CUDA kernel launch, aborting..."; \
    } \
    cuerr = cudaDeviceSynchronize(); \
    if(cuerr != cudaSuccess) { \
        printf("CUDA error in kernel execution: %s at kernel \"" #KernelCallInstruction "\"\n", cudaGetErrorString(cuerr)); \
            throw "error in CUDA kernel execution, aborting..."; \
    } \
}

__constant__ int d_cfI;
__constant__ int d_cfJ;
__constant__ int d_shapeYZ;
__constant__ int d_shapeZ;
__constant__ int d_gt_shapeYZ;
__constant__ int d_gt_shapeZ;
__constant__ int d_start_i;
__constant__ int d_start_j;
__constant__ int d_start_k;
__constant__ double d_tau;
__constant__ double d_h_x;
__constant__ double d_h_y;
__constant__ double d_h_z;
__constant__ double d_sqr_tau;

__constant__ double d_L_x;
__constant__ double d_L_y;
__constant__ double d_L_z;
__constant__ double d_a_t;

__device__
double laplacian(double *g, int index) {
    double center = g[index];
    return (g[index - d_cfI] - 2.0 * center + g[index + d_cfI]) / (d_h_x * d_h_x) +
           (g[index - d_cfJ] - 2.0 * center + g[index + d_cfJ]) / (d_h_y * d_h_y) +
           (g[index - 1] - 2.0 * center + g[index + 1]) / (d_h_z * d_h_z);
}

__global__
void cuda_step(double *grid, double *previous_1, double *previous_2) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i = idx / (d_shapeYZ) + 1;
    idx %= (d_shapeYZ);
    int j = idx / d_shapeZ + 1;
    int k = idx % d_shapeZ + 1;
    int index = i * d_cfI + j * d_cfJ + k;

    grid[index] = 2.0 * previous_1[index] - previous_2[index] +
                  d_sqr_tau * laplacian(previous_1, index);
}

__device__
double cuda_u(double t, double x, double y, double z) {
    return sin(2 * M_PI * x / d_L_x) * sin(M_PI * y / d_L_y) * sin(2 * M_PI * z / d_L_z) * cos(d_a_t * t);
}

__global__
void cuda_fillByGt(double *grid, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i = idx / (d_gt_shapeYZ);
    idx %= (d_gt_shapeYZ);
    int j = idx / d_gt_shapeZ;
    int k = idx % d_gt_shapeZ;
    int index = i * d_cfI + j * d_cfJ + k;
    grid[index] = cuda_u(
            d_tau * n,
            d_h_x * (d_start_i + i),
            d_h_y * (d_start_j + j),
            d_h_z * (d_start_k + k)
    );
}

void makeStepWithCuda(Grid3D &grid, Grid3D &previous_1, Grid3D &previous_2,
                      double h_x, double h_y, double h_z, double sqr_tau) {
    int blockSize = grid.shape[2] - 2;
    int gridInBlocks = (grid.shape[0] - 2) * (grid.shape[1] - 2);

    int shapeZ = grid.shape[2] - 2;
    int shapeYZ = (grid.shape[1] - 2) * shapeZ;

    size_t sizeInBytes = sizeof(double) * grid.size;

    double *d_grid;
    double *d_previous_1;
    double *d_previous_2;
    SAFE_CALL(cudaMalloc((void **) &d_grid, sizeInBytes));
    SAFE_CALL(cudaMalloc((void **) &d_previous_1, sizeInBytes));
    SAFE_CALL(cudaMalloc((void **) &d_previous_2, sizeInBytes));
    SAFE_CALL(cudaMemcpy(d_grid, grid.getFlatten().data(), sizeInBytes, cudaMemcpyHostToDevice));
    SAFE_CALL(cudaMemcpy(d_previous_1, previous_1.getFlatten().data(), sizeInBytes, cudaMemcpyHostToDevice));
    SAFE_CALL(cudaMemcpy(d_previous_2, previous_2.getFlatten().data(), sizeInBytes, cudaMemcpyHostToDevice));

    cudaMemcpyToSymbol(d_cfI, &grid._cfI, sizeof(int));
    cudaMemcpyToSymbol(d_cfJ, &grid._cfJ, sizeof(int));
    cudaMemcpyToSymbol(d_shapeYZ, &shapeYZ, sizeof(int));
    cudaMemcpyToSymbol(d_shapeZ, &shapeZ, sizeof(int));
    cudaMemcpyToSymbol(d_h_x, &h_x, sizeof(double));
    cudaMemcpyToSymbol(d_h_y, &h_y, sizeof(double));
    cudaMemcpyToSymbol(d_h_z, &h_z, sizeof(double));
    cudaMemcpyToSymbol(d_sqr_tau, &sqr_tau, sizeof(double));

    SAFE_KERNEL_CALL((cuda_step<<<gridInBlocks, blockSize>>>(d_grid, d_previous_1, d_previous_2)));

    SAFE_CALL(cudaMemcpy(grid.getFlatten().data(), d_grid, sizeInBytes, cudaMemcpyDeviceToHost));

    SAFE_CALL(cudaFree(d_grid));
    SAFE_CALL(cudaFree(d_previous_1));
    SAFE_CALL(cudaFree(d_previous_2));
}

void fillByGtWithCuda(Grid3D &grid, U u, int n, double tau, int start_i, int start_j, int start_k) {
    int blockSize = grid.shape[2];
    int gridInBlocks = grid.shape[0] * grid.shape[1];

    int shapeZ = grid.shape[2];
    int shapeYZ = grid.shape[1] * shapeZ;

    size_t sizeInBytes = sizeof(double) * grid.size;

    double *d_gt_grid;
    SAFE_CALL(cudaMalloc((void **) &d_gt_grid, sizeInBytes));

    cudaMemcpyToSymbol(d_cfI, &grid._cfI, sizeof(int));
    cudaMemcpyToSymbol(d_cfJ, &grid._cfJ, sizeof(int));
    cudaMemcpyToSymbol(d_gt_shapeYZ, &shapeYZ, sizeof(int));
    cudaMemcpyToSymbol(d_gt_shapeZ, &shapeZ, sizeof(int));
    cudaMemcpyToSymbol(d_start_i, &start_i, sizeof(int));
    cudaMemcpyToSymbol(d_start_j, &start_j, sizeof(int));
    cudaMemcpyToSymbol(d_start_k, &start_k, sizeof(int));
    cudaMemcpyToSymbol(d_L_x, &u.L_x, sizeof(double));
    cudaMemcpyToSymbol(d_L_y, &u.L_y, sizeof(double));
    cudaMemcpyToSymbol(d_L_z, &u.L_z, sizeof(double));
    cudaMemcpyToSymbol(d_a_t, &u.a_t, sizeof(double));
    cudaMemcpyToSymbol(d_tau, &tau, sizeof(double));

    SAFE_KERNEL_CALL((cuda_fillByGt<<<gridInBlocks, blockSize>>>(d_gt_grid, n)));

    SAFE_CALL(cudaMemcpy(grid.getFlatten().data(), d_gt_grid, sizeInBytes, cudaMemcpyDeviceToHost));

    SAFE_CALL(cudaFree(d_gt_grid));
}