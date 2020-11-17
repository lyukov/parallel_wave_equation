#include <cmath>
#include <fstream>
#include <iostream>

class Function4D {
   public:
    virtual double operator()(double t, double x, double y, double z) const = 0;
};

class Function3D {
   public:
    virtual double operator()(double x, double y, double z) const = 0;
};

class Grid {
    /* Raw data */
    double *raw;

   public:
    /* Grid size is K * N * N * N */
    const int N;
    const int K;

    /* Parameters of continuous problem */
    const double T;
    const double L_x;
    const double L_y;
    const double L_z;

    /* Grid steps */
    const double h_x;
    const double h_y;
    const double h_z;
    const double tau;

    const int size;

    Grid(double T, double L_x, double L_y, double L_z, int N, int K)
        : T(T),
          L_x(L_x),
          L_y(L_y),
          L_z(L_z),
          N(N),
          K(K),
          h_x(L_x / N),
          h_y(L_y / N),
          h_z(L_z / N),
          tau(T / K),
          size(K * N * N * N)  // May overflow ?
    {
        raw = new double[size];
    }

    ~Grid() {
        if (raw) delete raw;
    }

    double &operator()(int n, int i, int j, int k) {
        return raw[n * (N * N * N) + i * (N * N) + j * N + k];
    }

    void writeToFile(std::ofstream &outFile) const {
        outFile.write((char *)raw, size * sizeof(double));
    }
};

double laplasian(Grid g, int n, int i, int j, int k) {
    double center = g(n, i, j, k);
    return (g(n, i - 1, j, k) - 2.0 * center + g(n, i + 1, j, k)) / g.h_x +
           (g(n, i, j - 1, k) - 2.0 * center + g(n, i, j + 1, k)) / g.h_y +
           (g(n, i, j, k - 1) - 2.0 * center + g(n, i, j, k + 1)) / g.h_z;
}

class Solver {
    const Function4D *const u;
    const Function3D *const phi;

    const int N;
    const int K;

   public:
    Grid grid;

    Solver(double T, double L_x, double L_y, double L_z, int N, int K, Function4D *u, Function3D *phi)
        : N(N),
          K(K),
          u(u),
          phi(phi),
          grid(T, L_x, L_y, L_z, N, K) {}

    void init() {
        // Initialize zero level
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    grid(0, i, j, k) = (*phi)(
                        grid.h_x * i,
                        grid.h_y * j,
                        grid.h_z * k);
                }
            }
        }

        // Initialize first level
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    grid(1, i, j, k) = (*u)(  // TODO: initialize without using u
                        grid.tau,
                        grid.h_x * i,
                        grid.h_y * j,
                        grid.h_z * k);
                }
            }
        }
    }

    /**
     *  makeStep fills n-th layer of grid. It depends on two previous layers.
     * */
    void makeStep(int n) {
        if (n < 2 || n >= grid.T) {
            std::cerr << "Parameter n in makeStep must be between 2 and T. Actual value: " << n << std::endl;
        }

        // Inner nodes
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                for (int k = 1; k < N - 1; ++k) {
                    grid(n, i, j, k) = 2 * grid(n - 1, i, j, k) -
                                       grid(n - 2, i, j, k) +
                                       grid.tau * grid.tau * laplasian(grid, n - 1, i, j, k);
                }
            }
        }

        // Border nodes using first type condition
        // TODO: replace by my variant
        for (int i = 0; i < N; i += N - 1) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    grid(n, i, j, k) = 0;
                }
            }
        }
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; j += N - 1) {
                for (int k = 0; k < N; ++k) {
                    grid(n, i, j, k) = 0;
                }
            }
        }
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; k += N - 1) {
                    grid(n, i, j, k) = 0;
                }
            }
        }
    }

    void solve() {
        init();
        for (int step = 2; step < grid.T; ++step) {
            makeStep(step);
        }
    }
};

class U : Function4D {
    const double L_x;
    const double L_y;
    const double L_z;
    const double a_t;

   public:
    U(double L_x, double L_y, double L_z)
        : L_x(L_x),
          L_y(L_y),
          L_z(L_z),
          a_t(M_PI * sqrt(4.0 / (L_x * L_x) + 1.0 / (L_y * L_y) + 1.0 / (L_z * L_z))) {}

    double operator()(double t, double x, double y, double z) const {
        return sin(M_2_PI * x / L_x) * sin(M_PI * y / L_y) * sin(M_2_PI * z / L_z) * cos(a_t * t);
    }
};

class Phi : Function3D {
    double L_x;
    double L_y;
    double L_z;

   public:
    Phi(double L_x, double L_y, double L_z) {
        this->L_x = L_x;
        this->L_y = L_y;
        this->L_z = L_z;
    }

    double operator()(double x, double y, double z) const {
        return sin(M_2_PI * x / L_x) * sin(M_PI * y / L_y) * sin(M_2_PI * z / L_z);
    }
};

int main(int argc, char **argv) {
    if (argc <= 7) {
        std::cout << "Usage: prog L_x L_y L_z T N K out_file" << std::endl;
        return 0;
    }

    std::ofstream outFile(argv[7], std::ios::out | std::ios::binary);

    double L_x = atof(argv[1]);
    double L_y = atof(argv[2]);
    double L_z = atof(argv[3]);
    double T = atof(argv[4]);
    int N = atoi(argv[5]);
    int K = atoi(argv[6]);

    Phi phi(L_x, L_y, L_z);
    U u(L_x, L_y, L_z);

    Solver solver(K, L_x, L_y, L_z, N, K, (Function4D *)&u, (Function3D *)&phi);

    solver.solve();

    solver.grid.writeToFile(outFile);
    outFile.close();

    return 0;
}