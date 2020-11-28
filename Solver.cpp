#include "Solver.h"
#include "log.h"

Solver::Solver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi)
        : N(N),
          K(K),
          u(u),
          phi(phi),
          h_x(L_x / N),
          h_y(L_y / N),
          h_z(L_z / N),
          tau(T / K) {}

void Solver::init_0(Grid3D &grid, int start_i, int start_j, int start_k) const {
    // Initialize zero level
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
    LOG << "Level 0 initialized\n";
}

void Solver::init_1(Grid3D &grid, int start_i, int start_j, int start_k) const {
    // Initialize first level
    for (int i = 0; i < grid.shape[0]; ++i) {
        for (int j = 0; j < grid.shape[1]; ++j) {
            for (int k = 0; k < grid.shape[2]; ++k) {
                grid(i, j, k) = u(  // TODO: initialize without using u
                        tau,
                        h_x * (start_i + i),
                        h_y * (start_j + j),
                        h_z * (start_k + k)
                );
            }
        }
    }
    LOG << "Level 1 initialized\n";
}

double Solver::laplacian(const Grid3D &g, int i, int j, int k) const {
    double center = g(i, j, k);
    return (g(i - 1, j, k) - 2.0 * center + g(i + 1, j, k)) / h_x +
           (g(i, j - 1, k) - 2.0 * center + g(i, j + 1, k)) / h_y +
           (g(i, j, k - 1) - 2.0 * center + g(i, j, k + 1)) / h_z;
}

/** Fills n-th layer of grid. It depends on two previous layers. */
void Solver::fillInnerNodes(Grid3D &grid, const Grid3D &previous_1, const Grid3D &previous_2) const {
    // Inner nodes
    for (int i = 1; i < grid.shape[0] - 1; ++i) {
        for (int j = 1; j < grid.shape[1] - 1; ++j) {
            for (int k = 1; k < grid.shape[2] - 1; ++k) {
                grid(i, j, k) = 2 * previous_1(i, j, k) -
                                previous_2(i, j, k) +
                                 tau * tau * laplacian(previous_1, i, j, k);
            }
        }
    }
}

void Solver::fillByU(Grid3D &grid, int n) const {
    for (int i = 0; i < grid.shape[0]; ++i) {
        for (int j = 0; j < grid.shape[1]; ++j) {
            for (int k = 0; k < grid.shape[2]; ++k) {
                grid(i, j, k) = u(
                        tau * n,
                        h_x * i,
                        h_y * j,
                        h_z * k
                );
            }
        }
    }
}
