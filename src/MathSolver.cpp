#include "MathSolver.h"
#include "utils.h"

MathSolver::MathSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi)
        : N(N),
          K(K),
          u(u),
          phi(phi),
          h_x(L_x / N),
          h_y(L_y / N),
          h_z(L_z / N),
          tau(T / K) {}

void MathSolver::init_1(Grid3D &grid, int start_i, int start_j, int start_k) const {
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
    LOG_DEBUG << "Level 0 initialized\n";
}

void MathSolver::init_2(Grid3D &grid, int start_i, int start_j, int start_k) const {
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
    LOG_DEBUG << "Level 1 initialized\n";
}

double MathSolver::laplacian(const Grid3D &g, int i, int j, int k) const {
    double center = g(i, j, k);
    return (g(i - 1, j, k) - 2.0 * center + g(i + 1, j, k)) / h_x +
           (g(i, j - 1, k) - 2.0 * center + g(i, j + 1, k)) / h_y +
           (g(i, j, k - 1) - 2.0 * center + g(i, j, k + 1)) / h_z;
}

/** Fills n-th layer of grid. It depends on two previous layers. */
void MathSolver::makeStepForInnerNodes(Grid3D &grid, const Grid3D &previous_1, const Grid3D &previous_2) const {
    // Inner nodes
    for (int i = 1; i < grid.shape[0] - 1; ++i) {
        for (int j = 1; j < grid.shape[1] - 1; ++j) {
            for (int k = 1; k < grid.shape[2] - 1; ++k) {
                double result = 2 * previous_1(i, j, k) -
                                previous_2(i, j, k) +
                                tau * tau * laplacian(previous_1, i, j, k);
                if (abs(result) > 100) {
                    LOG << result << " " << i << " " << j << " " << k << endl;
                    LOG << previous_1(i, j, k) << " "
                        << previous_2(i, j, k) << " "
                        << laplacian(previous_1, i, j, k) << endl;
                }
                grid(i, j, k) = result;
            }
        }
    }
}

void MathSolver::fillByU(Grid3D &grid, int n, int start_i, int start_j, int start_k) const {
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

double MathSolver::C_norm_inner(const Grid3D &grid, const Grid3D &another) const {
    double c_norm = 0;
    for (int i = 1; i < grid.shape[0] - 1; ++i) {
        for (int j = 1; j < grid.shape[1] - 1; ++j) {
            for (int k = 1; k < grid.shape[2] - 1; ++k) {
                c_norm = max(
                        abs(grid(i, j, k) - another(i, j, k)),
                        c_norm
                );
            }
        }
    }
    return c_norm;
}

std::ostream &operator<<(std::ostream &out, const MathSolver &solver) {
    return out << "MathSolver: "
               << "N = " << solver.N << ", "
               << "K = " << solver.K << ", "
               << "h_x = " << solver.h_x << ", "
               << "h_y = " << solver.h_y << ", "
               << "h_z = " << solver.h_z << ", "
               << "tau = " << solver.tau;
}