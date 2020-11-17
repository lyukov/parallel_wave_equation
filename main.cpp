#include <cmath>
#include <iostream>

class Function4D {
   public:
    virtual double operator()(double t, double x, double y, double z) = 0;
};

class Function3D {
   public:
    virtual double operator()(double x, double y, double z) = 0;
};

class Grid {
    /* Raw data */
    double *raw;

   public:
    /* Grid size is K * N * N * N */
    int N;
    int K;

    /* Parameters of continuous problem */
    double T;
    double L_x;
    double L_y;
    double L_z;

    /* Grid steps */
    double h_x;
    double h_y;
    double h_z;
    double tau;

    Grid() {
        raw = 0;
    }

    Grid(double T, double L_x, double L_y, double L_z, int N, int K) {
        this->T = T;
        this->L_x = L_x;
        this->L_y = L_y;
        this->L_z = L_z;

        h_x = L_x / N;
        h_y = L_y / N;
        h_z = L_z / N;
        tau = T / K;

        long size = K * N * N * N;  // May overflow ?
        raw = new double[size];
    }

    ~Grid() {
        if (raw) delete raw;
    }

    double &operator()(int n, int i, int j, int k) {
        return raw[n * (N * N * N) + i * (N * N) + j * N + k];
    }

    double laplasian(int n, int i, int j, int k) {
        double center = (n, i, j, k);
        return ((n, i - 1, j, k) - 2.0 * center + (n, i + 1, j, k)) / h_x +
               ((n, i, j - 1, k) - 2.0 * center + (n, i, j + 1, k)) / h_y +
               ((n, i, j, k - 1) - 2.0 * center + (n, i, j, k + 1)) / h_z;
    }
};

class Solver {
    Function4D *u;
    Function3D *phi;

    int N = grid.N;
    int K;

   public:
    Grid grid;

    Solver(double T, double L_x, double L_y, double L_z, int N, int K, Function4D *u, Function3D *phi) {
        grid = Grid(T, L_x, L_y, L_z, N, K);
        this->u = u;
        this->phi = phi;
        this->N = N;
        this->K = K;
    }

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
     *  makeStep fills n-th layer of grid. It depends on two previous layers
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
                                       grid.tau * grid.tau * grid.laplasian(n - 1, i, j, k);
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
    double L_x;
    double L_y;
    double L_z;
    double a_t;

   public:
    U(double L_x, double L_y, double L_z) {
        this->L_x = L_x;
        this->L_y = L_y;
        this->L_z = L_z;
        a_t = M_PI * sqrt(4.0 / (L_x * L_x) + 1.0 / (L_y * L_y) + 1.0 / (L_z * L_z));
    }

    double operator()(double t, double x, double y, double z) {
        return sin(M_2_PI * x / L_x) * sin(M_PI * y / L_y) * sin(M_2_PI * z / L_z) * cos(a_t * t);
    }
};

class Phi : Function4D {
    double L_x;
    double L_y;
    double L_z;

   public:
    Phi(double L_x, double L_y, double L_z) {
        this->L_x = L_x;
        this->L_y = L_y;
        this->L_z = L_z;
    }

    double operator()(double x, double y, double z) {
        return sin(M_2_PI * x / L_x) * sin(M_PI * y / L_y) * sin(M_2_PI * z / L_z);
    }
};

int main(int argc, char **argv) {
    return 0;
}