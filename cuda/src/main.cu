#include <iostream>

using std::endl;

__global__ void sum_kernel(double *A, double *B, double *C, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) { return; }
    double a = A[idx];
    double b = B[idx];
    C[idx] = a + b;
}

int main(int argc, char **argv) {
    // Size of vectors
    int n = 100000;
    // Host vectors
    double *h_a, *h_b, *h_c;
    // Size, in bytes, of each vector
    size_t bytes = n * sizeof(double);
    // Allocate memory for each vector on host
    h_a = (double *) malloc(bytes);
    h_b = (double *) malloc(bytes);
    h_c = (double *) malloc(bytes);
    int i;
    // Initialize vectors on host
    for (i = 0; i < n; i++) {
        h_a[i] = sin(i) * sin(i);
        h_b[i] = cos(i) * cos(i);
    }

    // Device input vectors
    double *d_a, *d_b, *d_c;
    // Allocate memory for each vector on GPU
    cudaMalloc(&d_a, bytes);
    cudaMalloc(&d_b, bytes);
    cudaMalloc(&d_c, bytes);
    // Copy host vectors to device
    cudaMemcpy(d_a, h_a, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, h_b, bytes, cudaMemcpyHostToDevice);
    int blockSize = 1024;
    int gridSize = (n - 1) / blockSize + 1;
    // Execute the kernel
    sum_kernel<<<gridSize, blockSize>>>(d_a, d_b, d_c, n);
    // Copy array back to host
    cudaMemcpy(h_c, d_c, bytes, cudaMemcpyDeviceToHost);
    // Release device memory
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);

    double maxError = 0;
    for (int i = 0; i < n; ++i) {
        double error = abs(h_c[i] - 1.0);
        maxError = error > maxError ? error : maxError;
    }
    std::cout << "Max error = " << maxError << endl;
}