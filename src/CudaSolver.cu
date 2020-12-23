#include "CudaSolver.cuh"
#include <cstdlib>

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
__constant__ double d_h_x;
__constant__ double d_h_y;
__constant__ double d_h_z;
__constant__ double d_sqr_tau;

__device__
double &get_at(double *g, int i, int j, int k) {
    return g[i * d_cfI + j * d_cfJ + k];
}

__device__
double laplacian(double *g, int i, int j, int k) {
    double center = get_at(g, i, j, k);
    return (get_at(g, i - 1, j, k) - 2.0 * center + get_at(g, i + 1, j, k)) / (d_h_x * d_h_x) +
           (get_at(g, i, j - 1, k) - 2.0 * center + get_at(g, i, j + 1, k)) / (d_h_y * d_h_y) +
           (get_at(g, i, j, k - 1) - 2.0 * center + get_at(g, i, j, k + 1)) / (d_h_z * d_h_z);
}

__global__
void step(double *grid, double *previous_1, double *previous_2) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.z * blockDim.z + threadIdx.z + 1;
    get_at(grid, i, j, k) = 2.0 * get_at(previous_1, i, j, k) - get_at(previous_2, i, j, k) +
                            d_sqr_tau * laplacian(previous_1, i, j, k);
}

void makeStepWithCuda(Grid3D &grid, Grid3D &previous_1, Grid3D &previous_2,
                      double h_x, double h_y, double h_z, double sqr_tau) {
    dim3 blockSize = dim3(
            1,
            1,
            grid.shape[2] - 2
    );
    dim3 gridInBlocks = dim3(
            (grid.shape[0] - 2) / blockSize.x,
            (grid.shape[1] - 2) / blockSize.y,
            (grid.shape[2] - 2) / blockSize.z
    );

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
    cudaMemcpyToSymbol(d_h_x, &h_x, sizeof(double));
    cudaMemcpyToSymbol(d_h_y, &h_y, sizeof(double));
    cudaMemcpyToSymbol(d_h_z, &h_z, sizeof(double));
    cudaMemcpyToSymbol(d_sqr_tau, &sqr_tau, sizeof(double));

    step<<<gridInBlocks, blockSize>>>(d_grid, d_previous_1, d_previous_2);

    SAFE_CALL(cudaMemcpy(grid.getFlatten().data(), d_grid, sizeInBytes, cudaMemcpyDeviceToHost));
}