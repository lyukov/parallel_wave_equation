#include "CpuSolver.h"
#include "utils.h"
#include <omp.h>

CpuSolver::CpuSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi)
        : MathSolver(T, L_x, L_y, L_z, N, K, u, phi) {}

void CpuSolver::init_0(Grid3D &grid, int start_i, int start_j, int start_k) {
#pragma omp parallel for
    for (int i = 0; i < grid.shape[0]; ++i) {
        for (int j = 0; j < grid.shape[1]; ++j) {
            for (int k = 0; k < grid.shape[2]; ++k) {
                grid(i, j, k) = phi(
                        h_x * (start_i + i),
                        h_y * (start_j + j),
                        h_z * (start_k + k)
                );
            }
        }
    }
}

void CpuSolver::init_1(Grid3D &grid, Grid3D &previous) {
#pragma omp parallel for
    for (int i = 1; i < grid.shape[0] - 1; ++i) {
        for (int j = 1; j < grid.shape[1] - 1; ++j) {
            for (int k = 1; k < grid.shape[2] - 1; ++k) {
                grid(i, j, k) = previous(i, j, k) +
                                0.5 * sqr(tau) * laplacian(previous, i, j, k);
            }
        }
    }
}

double CpuSolver::laplacian(const Grid3D &g, int i, int j, int k) const {
    double center = g(i, j, k);
    return (g(i - 1, j, k) - 2.0 * center + g(i + 1, j, k)) / sqr(h_x) +
           (g(i, j - 1, k) - 2.0 * center + g(i, j + 1, k)) / sqr(h_y) +
           (g(i, j, k - 1) - 2.0 * center + g(i, j, k + 1)) / sqr(h_z);
}

/** Fills n-th layer of grid. It depends on two previous layers. */
void CpuSolver::makeStepForInnerNodes(Grid3D &grid, Grid3D &previous_1, Grid3D &previous_2) {
#pragma omp parallel for
    for (int i = 1; i < grid.shape[0] - 1; ++i) {
        for (int j = 1; j < grid.shape[1] - 1; ++j) {
            for (int k = 1; k < grid.shape[2] - 1; ++k) {
                double result = 2.0 * previous_1(i, j, k) - previous_2(i, j, k) +
                                sqr(tau) * laplacian(previous_1, i, j, k);
                if (grid(i, j, k) != result) {
                    LOG << "unmatch: " << grid(i, j, k) << " " << result << endl;
                }
            }
        }
    }
}

void CpuSolver::fillByGroundTruth(Grid3D &grid, int n, int start_i, int start_j, int start_k) {
#pragma omp parallel for
    for (int i = 0; i < grid.shape[0]; ++i) {
        for (int j = 0; j < grid.shape[1]; ++j) {
            for (int k = 0; k < grid.shape[2]; ++k) {
                grid(i, j, k) = u(
                        tau * n,
                        h_x * (start_i + i),
                        h_y * (start_j + j),
                        h_z * (start_k + k)
                );
            }
        }
    }
}

double CpuSolver::maxAbsoluteErrorInner(Grid3D &grid, Grid3D &another) {
    double c_norm = 0;
    for (int i = 1; i < grid.shape[0] - 1; ++i) {
        for (int j = 1; j < grid.shape[1] - 1; ++j) {
            for (int k = 1; k < grid.shape[2] - 1; ++k) {
                c_norm = max(
                        std::abs(grid(i, j, k) - another(i, j, k)),
                        c_norm
                );
            }
        }
    }
    return c_norm;
}

double CpuSolver::sumSquaredErrorInner(Grid3D &grid, Grid3D &another) {
    double sum = 0;
    for (int i = 1; i < grid.shape[0] - 1; ++i) {
        for (int j = 1; j < grid.shape[1] - 1; ++j) {
            for (int k = 1; k < grid.shape[2] - 1; ++k) {
                sum += sqr(grid(i, j, k) - another(i, j, k));
            }
        }
    }
    return sum;
}
