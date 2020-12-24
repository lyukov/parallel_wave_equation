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
} blockIdx, blockDim, threadIdx;

struct cudaError_t {
    bool operator!=(cudaError_t &other);
} cudaSuccess;

template<typename T>
cudaError_t cudaMemcpyToSymbol(T &symbol, const void *src, size_t count);

cudaError_t cudaMemcpy(void *dst, const void *src, size_t count, enum cudaMemcpyKind kind);

cudaError_t cudaMemset(void *dst, int value, size_t count);

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
    return (g[index - d_shapeYZ] - 2.0 * center + g[index + d_shapeYZ]) / (d_h_x * d_h_x) +
           (g[index - d_shapeZ] - 2.0 * center + g[index + d_shapeZ]) / (d_h_y * d_h_y) +
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
    return i * d_shapeYZ + j * d_shapeZ + k;
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

__global__
void cuda_c1(double *left, double *right, double *result) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    coords_inner(idx, i, j, k);
    int index = flat_index(i, j, k);
    result[index] = abs(left[index] - right[index]);
}

__global__
void cuda_squared_error(double *left, double *right, double *result) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i, j, k;
    coords_inner(idx, i, j, k);
    int index = flat_index(i, j, k);
    result[index] = pow(left[index] - right[index], 2);
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

__global__
void cuda_get_slice(int c0, int c1, int c2, double *slice, double *grid) {
    int grid_idx = c0 + blockIdx.x * c1 + threadIdx.x * c2;
    int slice_idx = blockIdx.x * blockDim.x + threadIdx.x;
    slice[slice_idx] = grid[grid_idx];
}

__global__
void cuda_set_slice(int c0, int c1, int c2, double *slice, double *grid) {
    int grid_idx = c0 + blockIdx.x * c1 + threadIdx.x * c2;
    int slice_idx = blockIdx.x * blockDim.x + threadIdx.x;
    grid[grid_idx] = slice[slice_idx];
}

CudaSolver::CudaSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi,
                       int shapeX, int shapeY, int shapeZ)
        : MathSolver(T, L_x, L_y, L_z, N, K, u, phi),
          sizeInBytes(sizeof(double) * shapeX * shapeY * shapeZ),
          flatSize(shapeX * shapeY * shapeZ),
          gridSizeFull(shapeX * shapeY),
          blockSizeFull(shapeZ),
          gridSizeInner((shapeX - 2) * (shapeY - 2)),
          blockSizeInner(shapeZ - 2),
          grid3D(shapeX, shapeY, shapeZ) {

    cudaMemcpyToSymbol(d_h_x, &h_x, sizeof(double));
    cudaMemcpyToSymbol(d_h_y, &h_y, sizeof(double));
    cudaMemcpyToSymbol(d_h_z, &h_z, sizeof(double));
    cudaMemcpyToSymbol(d_tau, &tau, sizeof(double));
    cudaMemcpyToSymbol(d_L_x, &u.L_x, sizeof(double));
    cudaMemcpyToSymbol(d_L_y, &u.L_y, sizeof(double));
    cudaMemcpyToSymbol(d_L_z, &u.L_z, sizeof(double));
    cudaMemcpyToSymbol(d_a_t, &u.a_t, sizeof(double));

    int shapeYZ = shapeY * shapeZ;
    int shapeZ_inner = shapeZ - 2;
    int shapeYZ_inner = (shapeY - 2) * shapeZ_inner;
    cudaMemcpyToSymbol(d_shapeYZ, &shapeYZ, sizeof(int));
    cudaMemcpyToSymbol(d_shapeZ, &shapeZ, sizeof(int));
    cudaMemcpyToSymbol(d_shapeYZ_inner, &shapeYZ_inner, sizeof(int));
    cudaMemcpyToSymbol(d_shapeZ_inner, &shapeZ_inner, sizeof(int));

    d_grids.resize(N_GRIDS);
    for (int i = 0; i < N_GRIDS; ++i) {
        SAFE_CALL(cudaMalloc((void **) &d_grids[i], sizeInBytes));
    }
    SAFE_CALL(cudaMalloc((void **) &d_groundTruth, sizeInBytes));
    SAFE_CALL(cudaMalloc((void **) &d_errorC1, sizeInBytes));
    SAFE_CALL(cudaMalloc((void **) &d_errorMSE, sizeInBytes));
    SAFE_CALL(cudaMemset(d_errorC1, 0, sizeInBytes));
    SAFE_CALL(cudaMemset(d_errorMSE, 0, sizeInBytes));

    int maxSliceSize = max(grid3D.getSliseSize(0), max(grid3D.getSliseSize(1), grid3D.getSliseSize(2)));
    h_slice.resize(maxSliceSize);
    SAFE_CALL(cudaMalloc((void **) &d_slice, sizeof(double) * maxSliceSize));
}

CudaSolver::~CudaSolver() {
    for (int i = 0; i < N_GRIDS; ++i) {
        SAFE_CALL(cudaFree(d_grids[i]));
    }
    SAFE_CALL(cudaFree(d_groundTruth));
    SAFE_CALL(cudaFree(d_errorC1));
    SAFE_CALL(cudaFree(d_errorMSE));
    SAFE_CALL(cudaFree(d_slice));
}

void CudaSolver::init_0(int start_i, int start_j, int start_k) {
    double *d_grid = d_grids[0];
    cudaMemcpyToSymbol(d_start_i, &start_i, sizeof(int));
    cudaMemcpyToSymbol(d_start_j, &start_j, sizeof(int));
    cudaMemcpyToSymbol(d_start_k, &start_k, sizeof(int));
    SAFE_KERNEL_CALL((cuda_init0<<<gridSizeInner, blockSizeInner>>>(d_grid)));
}

void CudaSolver::init_1() {
    double *d_grid = d_grids[1];
    double *d_previous = d_grids[0];
    SAFE_KERNEL_CALL((cuda_init1<<<gridSizeInner, blockSizeInner>>>(d_grid, d_previous)));
}

void CudaSolver::makeStepForInnerNodes(int n) {
    double *d_grid = d_grids[n % N_GRIDS];
    double *d_previous_1 = d_grids[(n - 1) % N_GRIDS];
    double *d_previous_2 = d_grids[(n - 2) % N_GRIDS];
    SAFE_KERNEL_CALL((cuda_step<<<gridSizeInner, blockSizeInner>>>(d_grid, d_previous_1, d_previous_2)));
}

void CudaSolver::updateGroundTruth(int n, int start_i, int start_j, int start_k) {
    cudaMemcpyToSymbol(d_start_i, &start_i, sizeof(int));
    cudaMemcpyToSymbol(d_start_j, &start_j, sizeof(int));
    cudaMemcpyToSymbol(d_start_k, &start_k, sizeof(int));
    SAFE_KERNEL_CALL((cuda_fillByGt<<<gridSizeFull, blockSizeFull>>>(d_groundTruth, n)));
}

double CudaSolver::maxAbsoluteErrorInner(int n) {
    double *d_grid = d_grids[n % N_GRIDS];
    SAFE_KERNEL_CALL((cuda_c1<<<gridSizeInner, blockSizeInner>>>(d_grid, d_groundTruth, d_errorC1)));
    return thrust::reduce(thrust::device, d_errorC1, d_errorC1 + flatSize, 0.0, thrust::maximum<double>());
}

double CudaSolver::sumSquaredErrorInner(int n) {
    double *d_grid = d_grids[n % N_GRIDS];
    SAFE_KERNEL_CALL((cuda_squared_error<<<gridSizeInner, blockSizeInner>>>(d_grid, d_groundTruth, d_errorMSE)));
    return thrust::reduce(thrust::device, d_errorMSE, d_errorMSE + flatSize, 0.0, thrust::plus<double>());
}

void CudaSolver::getSliceParams(int axis, int &c0, int &c1, int &c2, int &gridSize, int &blockSize) const {
    int cf[3] = {grid3D.shape[1] * grid3D.shape[2], grid3D.shape[2], 1};
    c0 = cf[axis % 3];
    c1 = cf[(axis + 1) % 3];
    c2 = cf[(axis + 2) % 3];
    gridSize = grid3D.shape[(axis + 1) % N_GRIDS];
    blockSize = grid3D.shape[(axis + 2) % N_GRIDS];
}

int CudaSolver::getSliceSize(int axis) {
    return grid3D.getSliceSize(axis);
}

std::vector<double> CudaSolver::getSlice(int n, int index, int axis) {
    int c0, c1, c2, gridSize, blockSize;
    getSliceParams(axis, c0, c1, c2, gridSize, blockSize);
    SAFE_KERNEL_CALL((cuda_get_slice<<<gridSize, blockSize>>>(
            c0 * index, c1, c2, d_slice, getCurrentState(n)
    )));
    SAFE_CALL(cudaMemcpy(h_slice.data(), d_slice, getSliceSize(axis) * sizeof(double), cudaMemcpyDeviceToHost));
    return std::vector<double>(h_slice.begin(), h_slice.begin() + getSliceSize(axis));
}

void CudaSolver::setSlice(int n, int index, int axis, std::vector<double> &slice) {
    int c0, c1, c2, gridSize, blockSize;
    getSliceParams(axis, c0, c1, c2, gridSize, blockSize);
    SAFE_CALL(cudaMemcpy(d_slice, slice.data(), getSliceSize(axis) * sizeof(double), cudaMemcpyHostToDevice));
    SAFE_KERNEL_CALL((cuda_set_slice<<<gridSize, blockSize>>>(
            c0 * index, c1, c2, d_slice, getCurrentState(n)
    )));
}

void CudaSolver::setZeros(int n, int index, int axis) {
    int c0, c1, c2, gridSize, blockSize;
    getSliceParams(axis, c0, c1, c2, gridSize, blockSize);
    SAFE_CALL(cudaMemset(d_slice, 0, getSliceSize(axis) * sizeof(double)));
    SAFE_KERNEL_CALL((cuda_set_slice<<<gridSize, blockSize>>>(
            c0 * index, c1, c2, d_slice, getCurrentState(n)
    )));
}

double *CudaSolver::getCurrentState(int n) {
    return d_grids[n % N_GRIDS];
}
