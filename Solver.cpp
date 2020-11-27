#include "Solver.h"
#include "log.h"

Solver::Solver(double T, double L_x, double L_y, double L_z,
               int N, int K, Function4D *u, Function3D *phi)
        : N(N),
          K(K),
          u(u),
          phi(phi),
          T(T),
          L_x(L_x),
          L_y(L_y),
          L_z(L_z),
          h_x(L_x / N),
          h_y(L_y / N),
          h_z(L_z / N),
          tau(T / K) {}

void Solver::init_0(Grid3D &grid) {
    // Initialize zero level
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            for (int k = 0; k <= N; ++k) {
                grid(i, j, k) = (*phi)(
                        h_x * i,
                        h_y * j,
                        h_z * k);
            }
        }
    }
    LOG << "Level 0 initialized\n";
}

void Solver::init_1(Grid3D &grid) {
    // Initialize first level
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            for (int k = 0; k <= N; ++k) {
                grid(i, j, k) = (*u)(  // TODO: initialize without using u
                        tau,
                        h_x * i,
                        h_y * j,
                        h_z * k);
            }
        }
    }
    LOG << "Level 1 initialized\n";
}

double Solver::laplasian(const Grid3D &g, int i, int j, int k) const {
    double center = g(i, j, k);
    return (g(i - 1, j, k) - 2.0 * center + g(i + 1, j, k)) / h_x +
           (g(i, j - 1, k) - 2.0 * center + g(i, j + 1, k)) / h_y +
           (g(i, j, k - 1) - 2.0 * center + g(i, j, k + 1)) / h_z;
}

/** makeStep fills n-th layer of grid. It depends on two previous layers. */
void Solver::makeStep(const int n,
                      Grid3D &layer,
                      const Grid3D &previous_1,
                      const Grid3D &previous_2) const {
    if (n < 2 || n >= K) {
        LOG_ERR << "Parameter n in makeStep must be between 2 and T. Actual value: " << n << std::endl;
    }

    // Inner nodes
    for (int i = 1; i <= N - 1; ++i) {
        for (int j = 1; j <= N - 1; ++j) {
            for (int k = 1; k <= N - 1; ++k) {
                layer(i, j, k) = 2 * previous_1(i, j, k) -
                                 previous_2(i, j, k) +
                                 tau * tau * laplasian(previous_1, i, j, k);
            }
        }
    }

    // Border nodes using first type condition
    // TODO: replace by my variant
    for (int i = 0; i <= N; i += N) {
        for (int j = 0; j <= N; ++j) {
            for (int k = 0; k <= N; ++k) {
                layer(i, j, k) = 0;
            }
        }
    }
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; j += N) {
            for (int k = 0; k <= N; ++k) {
                layer(i, j, k) = 0;
            }
        }
    }
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            for (int k = 0; k <= N; k += N) {
                layer(i, j, k) = 0;
            }
        }
    }

    if ((n + 1) % 50 == 0) {
        LOG << "Step " << n + 1 << " completed" << std::endl;
    }
}

Grid3D *Solver::solve() {
    LOG << "Start solve()\n";
    Grid3D *grids[3] = {new Grid3D(N), new Grid3D(N), new Grid3D(N)};
    LOG << "Grids created\n";
    init_0(*grids[0]);
    init_1(*grids[1]);
    for (int step = 2; step <= K; ++step) {
        makeStep(
                step,
                *grids[step % 3],
                *grids[(step - 1) % 3],
                *grids[(step - 2) % 3]);
    }
    delete grids[(K - 1) % 3];
    delete grids[(K - 2) % 3];
    return grids[K % 3];
}

void Solver::fillByF(Grid3D &grid, Function4D *f, int n) {
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            for (int k = 0; k <= N; ++k) {
                grid(i, j, k) = (*f)(
                        tau * n,
                        h_x * i,
                        h_y * j,
                        h_z * k);
            }
        }
    }
}
