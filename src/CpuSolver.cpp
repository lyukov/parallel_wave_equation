#include "CpuSolver.h"
#include "utils.h"
#include <omp.h>

CpuSolver::CpuSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi,
                     int shapeX, int shapeY, int shapeZ)
        : MathSolver(T, L_x, L_y, L_z, N, K, u, phi),
          groundTruth(shapeX, shapeY, shapeZ) {
    for (int i = 0; i < N_GRIDS; ++i) {
        grids.push_back(
                Grid3D(shapeX, shapeY, shapeZ)
        );
    }
}

void CpuSolver::init_0(int start_i, int start_j, int start_k) {
    Grid3D &grid = grids[0];
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

void CpuSolver::init_1() {
    Grid3D &grid = grids[1];
    Grid3D &previous = grids[0];
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
void CpuSolver::makeStepForInnerNodes(int n) {
    Grid3D &grid = grids[n % N_GRIDS];
    Grid3D &previous_1 = grids[(n - 1) % N_GRIDS];
    Grid3D &previous_2 = grids[(n - 2) % N_GRIDS];
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

void CpuSolver::updateGroundTruth(int n, int start_i, int start_j, int start_k) {
#pragma omp parallel for
    for (int i = 0; i < groundTruth.shape[0]; ++i) {
        for (int j = 0; j < groundTruth.shape[1]; ++j) {
            for (int k = 0; k < groundTruth.shape[2]; ++k) {
                groundTruth(i, j, k) = u(
                        tau * n,
                        h_x * (start_i + i),
                        h_y * (start_j + j),
                        h_z * (start_k + k)
                );
            }
        }
    }
}

double CpuSolver::maxAbsoluteErrorInner(int n) {
    Grid3D &grid = grids[n % N_GRIDS];
    double c_norm = 0;
#pragma omp parallel for reduction(max : c_norm)
    for (int i = 1; i < grid.shape[0] - 1; ++i) {
        for (int j = 1; j < grid.shape[1] - 1; ++j) {
            for (int k = 1; k < grid.shape[2] - 1; ++k) {
                c_norm = std::max(
                        std::abs(grid(i, j, k) - groundTruth(i, j, k)),
                        c_norm
                );
            }
        }
    }
    return c_norm;
}

double CpuSolver::sumSquaredErrorInner(int n) {
    Grid3D &grid = grids[n % N_GRIDS];
    double sum = 0;
#pragma omp parallel for reduction(+ : sum)
    for (int i = 1; i < grid.shape[0] - 1; ++i) {
        for (int j = 1; j < grid.shape[1] - 1; ++j) {
            for (int k = 1; k < grid.shape[2] - 1; ++k) {
                sum += sqr(grid(i, j, k) - groundTruth(i, j, k));
            }
        }
    }
    return sum;
}

Grid3D &CpuSolver::getCurrentState(int n) {
    return grids[n % N_GRIDS];
}

std::vector<double> CpuSolver::getSlice(int n, int index, int axis) {
    return getCurrentState(n).getSlice(index, axis);
}

int CpuSolver::getSliceSize(int axis) {
    return groundTruth.getSliceSize(axis);
}

void CpuSolver::setSlice(int n, int index, int axis, std::vector<double> &slice) {
    getCurrentState(n).setSlice(index, axis, slice);
}

void CpuSolver::setZeros(int n, int index, int axis) {
    getCurrentState(n).setZeros(index, axis);
}

double CpuSolver::maxGroundTruth() {
    double res = 0.0;
#pragma omp parallel for reduction(max : res)
    for (int i = 0; i < groundTruth.shape[0]; ++i) {
        for (int j = 0; j < groundTruth.shape[1]; ++j) {
            for (int k = 0; k < groundTruth.shape[2]; ++k) {
                res = std::max(res, std::abs(groundTruth(i, j, k)));
            }
        }
    }
    return res;
}
