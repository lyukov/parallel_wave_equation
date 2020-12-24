#include "CudaSolver.cuh"
#include "utils.h"
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <cstdlib>
#include <cuda.h>
#include <cmath>

#ifndef __NVCC__

#define __constant__
#define __device__
#define __global__
#define __host__
typedef unsigned long size_t;
enum cudaMemcpyKind {
    cudaMemcpyHostToHost,
    cudaMemcpyHostToDevice,
    cudaMemcpyDeviceToHost,
    cudaMemcpyDeviceToDevice
};

struct dim3 {
    int x, y, z;

    dim3(int x = 1, int y = 1, int z = 1);
} blockIdx, blockDim, threadIdx; //

struct cudaError_t {
    bool operator!=(cudaError_t &other);
} cudaSuccess;

template<typename T>
cudaError_t cudaMemcpyToSymbol(T &symbol, const void *src, size_t count);

cudaError_t cudaMemcpy(void *dst, const void *src, size_t count, enum cudaMemcpyKind kind);

char *cudaGetErrorString(cudaError_t err);

int printf(const char *format, ...);

cudaError_t cudaFree(void *p);

cudaError_t cudaMalloc(void **dst, size_t size);

cudaError_t cudaGetLastError();

cudaError_t cudaDeviceSynchronize();

typedef int cudaStream_t;

cudaError_t *cudaConfigureCall(dim3 gridDim, dim3 blockDim, size_t sharedMem = 0, cudaStream_t stream = 0);

#endif

#define SAFE_CALL(Call) { \
    cudaError_t cuerr = Call; \
    if(cuerr != cudaSuccess) { \
        printf("CUDA error: %s at call \"" #Call "\"\n", cudaGetErrorString(cuerr)); \
            throw "error in CUDA API function, aborting..."; \
    } \
}

#define SAFE_KERNEL_CALL(KernelCall) { \
    KernelCall; \
    cudaError_t cuerr = cudaGetLastError(); \
    if(cuerr != cudaSuccess) { \
        printf("CUDA error in kernel launch: %s at kernel \"" #KernelCall "\"\n", cudaGetErrorString(cuerr)); \
            throw "error in CUDA kernel launch, aborting..."; \
    } \
    cuerr = cudaDeviceSynchronize(); \
    if(cuerr != cudaSuccess) { \
        printf("CUDA error in kernel execution: %s at kernel \"" #KernelCall "\"\n", cudaGetErrorString(cuerr)); \
            throw "error in CUDA kernel execution, aborting..."; \
    } \
}

__constant__ int d_cfI;
__constant__ int d_cfJ;
__constant__ int d_shapeYZ;
__constant__ int d_shapeZ;
__constant__ int d_shapeYZ_inner;
__constant__ int d_shapeZ_inner;
__constant__ int d_start_i;
__constant__ int d_start_j;
__constant__ int d_start_k;
__constant__ double d_tau;
__constant__ double d_h_x;
__constant__ double d_h_y;
__constant__ double d_h_z;

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

__device__
void coords_inner(int idx, int &i, int &j, int &k) {
    i = idx / (d_shapeYZ_inner) + 1;
    idx %= (d_shapeYZ_inner);
    j = idx / d_shapeZ_inner + 1;
    k = idx % d_shapeZ_inner + 1;
}

__device__
void coords_full(int idx, int &i, int &j, int &k) {
    i = idx / (d_shapeYZ);
    idx %= (d_shapeYZ);
    j = idx / d_shapeZ;
    k = idx % d_shapeZ;
}

__device__
int flat_index(int i, int j, int k) {
    return i * d_cfI + j * d_cfJ + k;
}

__global__
void cuda_step(double *grid, double *previous_1, double *previous_2) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    coords_inner(idx, i, j, k);
    int index = flat_index(i, j, k);

    grid[index] = 2.0 * previous_1[index] - previous_2[index] +
                  d_tau * d_tau * laplacian(previous_1, index);
}

__device__
double cuda_u(double t, double x, double y, double z) {
    return sin(2 * M_PI * x / d_L_x) * sin(M_PI * y / d_L_y) * sin(2 * M_PI * z / d_L_z) * cos(d_a_t * t);
}

__device__
double cuda_phi(double x, double y, double z) {
    return sin(2 * M_PI * x / d_L_x) * sin(M_PI * y / d_L_y) * sin(2 * M_PI * z / d_L_z);
}

__global__
void cuda_fillByGt(double *grid, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    coords_full(idx, i, j, k);
    int index = flat_index(i, j, k);
    grid[index] = cuda_u(
            d_tau * n,
            d_h_x * (d_start_i + i),
            d_h_y * (d_start_j + j),
            d_h_z * (d_start_k + k)
    );
}

__global__
void cuda_init0(double *grid) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    coords_inner(idx, i, j, k);
    int index = flat_index(i, j, k);
    grid[index] = cuda_phi(
            d_h_x * (d_start_i + i),
            d_h_y * (d_start_j + j),
            d_h_z * (d_start_k + k)
    );
}

__global__
void cuda_init1(double *grid, double *previous) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    coords_inner(idx, i, j, k);
    int index = flat_index(i, j, k);
    grid[index] = previous[index] + 0.5 * d_tau * d_tau * laplacian(previous, index);
}

CudaSolver::CudaSolver(
        double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi, Grid3D &grid
) : CpuSolver(T, L_x, L_y, L_z, N, K, u, phi),
    sizeInBytes(sizeof(double) * grid.size),
    flatSize(grid.size),
    gridSizeFull(grid.shape[0] * grid.shape[1]),
    blockSizeFull(grid.shape[2]),
    gridSizeInner((grid.shape[0] - 2) * (grid.shape[1] - 2)),
    blockSizeInner(grid.shape[2] - 2) {

    cudaMemcpyToSymbol(d_h_x, &h_x, sizeof(double));
    cudaMemcpyToSymbol(d_h_y, &h_y, sizeof(double));
    cudaMemcpyToSymbol(d_h_z, &h_z, sizeof(double));
    cudaMemcpyToSymbol(d_tau, &tau, sizeof(double));
    cudaMemcpyToSymbol(d_L_x, &u.L_x, sizeof(double));
    cudaMemcpyToSymbol(d_L_y, &u.L_y, sizeof(double));
    cudaMemcpyToSymbol(d_L_z, &u.L_z, sizeof(double));
    cudaMemcpyToSymbol(d_a_t, &u.a_t, sizeof(double));

    cudaMemcpyToSymbol(d_cfI, &grid._cfI, sizeof(int));
    cudaMemcpyToSymbol(d_cfJ, &grid._cfJ, sizeof(int));
    int shapeZ = grid.shape[2];
    int shapeYZ = grid.shape[1] * shapeZ;
    int shapeZ_inner = grid.shape[2] - 2;
    int shapeYZ_inner = (grid.shape[1] - 2) * shapeZ_inner;
    cudaMemcpyToSymbol(d_shapeYZ, &shapeYZ, sizeof(int));
    cudaMemcpyToSymbol(d_shapeZ, &shapeZ, sizeof(int));
    cudaMemcpyToSymbol(d_shapeYZ_inner, &shapeYZ_inner, sizeof(int));
    cudaMemcpyToSymbol(d_shapeZ_inner, &shapeZ_inner, sizeof(int));
}

void CudaSolver::init_0(Grid3D &grid, int start_i, int start_j, int start_k) {
    double *d_gt_grid;
    SAFE_CALL(cudaMalloc((void **) &d_gt_grid, sizeInBytes));
    cudaMemcpyToSymbol(d_start_i, &start_i, sizeof(int));
    cudaMemcpyToSymbol(d_start_j, &start_j, sizeof(int));
    cudaMemcpyToSymbol(d_start_k, &start_k, sizeof(int));
    SAFE_KERNEL_CALL((cuda_init0<<<gridSizeInner, blockSizeInner>>>(d_gt_grid)));
    SAFE_CALL(cudaMemcpy(grid.getFlatten().data(), d_gt_grid, sizeInBytes, cudaMemcpyDeviceToHost));
    SAFE_CALL(cudaFree(d_gt_grid));
}

void CudaSolver::init_1(Grid3D &grid, Grid3D &previous) {
//    CpuSolver::init_1(grid, previous);
    double *d_gt_grid;
    double *d_previous;
    SAFE_CALL(cudaMalloc((void **) &d_gt_grid, sizeInBytes));
    SAFE_CALL(cudaMalloc((void **) &d_previous, sizeInBytes));
    SAFE_CALL(cudaMemcpy(d_previous, previous.getFlatten().data(), sizeInBytes, cudaMemcpyHostToDevice));
    SAFE_KERNEL_CALL((cuda_init1<<<gridSizeInner, blockSizeInner>>>(d_gt_grid, d_previous)));
    SAFE_CALL(cudaMemcpy(grid.getFlatten().data(), d_gt_grid, sizeInBytes, cudaMemcpyDeviceToHost));
    SAFE_CALL(cudaFree(d_gt_grid));
    SAFE_CALL(cudaFree(d_previous));
}

void CudaSolver::makeStepForInnerNodes(Grid3D &grid, Grid3D &previous_1, Grid3D &previous_2) {
    double *d_grid;
    double *d_previous_1;
    double *d_previous_2;
    SAFE_CALL(cudaMalloc((void **) &d_grid, sizeInBytes));
    SAFE_CALL(cudaMalloc((void **) &d_previous_1, sizeInBytes));
    SAFE_CALL(cudaMalloc((void **) &d_previous_2, sizeInBytes));
    SAFE_CALL(cudaMemcpy(d_previous_1, previous_1.getFlatten().data(), sizeInBytes, cudaMemcpyHostToDevice));
    SAFE_CALL(cudaMemcpy(d_previous_2, previous_2.getFlatten().data(), sizeInBytes, cudaMemcpyHostToDevice));
    SAFE_KERNEL_CALL((cuda_step<<<gridSizeInner, blockSizeInner>>>(d_grid, d_previous_1, d_previous_2)));
    SAFE_CALL(cudaMemcpy(grid.getFlatten().data(), d_grid, sizeInBytes, cudaMemcpyDeviceToHost));
    SAFE_CALL(cudaFree(d_grid));
    SAFE_CALL(cudaFree(d_previous_1));
    SAFE_CALL(cudaFree(d_previous_2));
}

void CudaSolver::fillByGroundTruth(Grid3D &grid, int n, int start_i, int start_j, int start_k) {
    double *d_gt_grid;
    SAFE_CALL(cudaMalloc((void **) &d_gt_grid, sizeInBytes));
    cudaMemcpyToSymbol(d_start_i, &start_i, sizeof(int));
    cudaMemcpyToSymbol(d_start_j, &start_j, sizeof(int));
    cudaMemcpyToSymbol(d_start_k, &start_k, sizeof(int));
    SAFE_KERNEL_CALL((cuda_fillByGt<<<gridSizeFull, blockSizeFull>>>(d_gt_grid, n)));
    SAFE_CALL(cudaMemcpy(grid.getFlatten().data(), d_gt_grid, sizeInBytes, cudaMemcpyDeviceToHost));
    SAFE_CALL(cudaFree(d_gt_grid));
}

struct cuda_c1 {
    __host__ __device__
    double operator()(const double &x1, const double &x2) const {
        return abs(x1 - x2);
    }
};

double CudaSolver::maxAbsoluteErrorInner(Grid3D &grid, Grid3D &another) {
    //return CpuSolver::maxAbsoluteErrorInner(grid, another);
    double *d_grid, *d_another, *d_error;
    SAFE_CALL(cudaMalloc((void **) &d_grid, sizeInBytes));
    SAFE_CALL(cudaMalloc((void **) &d_another, sizeInBytes));
    SAFE_CALL(cudaMalloc((void **) &d_error, sizeInBytes));
    SAFE_CALL(cudaMemcpy(d_grid, grid.getFlatten().data(), sizeInBytes, cudaMemcpyHostToDevice));
    SAFE_CALL(cudaMemcpy(d_another, another.getFlatten().data(), sizeInBytes, cudaMemcpyHostToDevice));
    SAFE_KERNEL_CALL(thrust::transform(d_grid, d_grid + flatSize, d_another, d_error, cuda_c1()));
    LOG << "Transform finished" << endl;
    double error = thrust::reduce(thrust::device, d_error, d_error + flatSize, 0.0, thrust::maximum<double>());
    LOG << "Reduce finished" << endl;
    SAFE_CALL(cudaFree(d_grid));
    SAFE_CALL(cudaFree(d_another));
    SAFE_CALL(cudaFree(d_error));
    return error;
}
