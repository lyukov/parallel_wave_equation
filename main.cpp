#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>

using std::endl;

std::string getTimestamp() {
    std::time_t time = std::time(nullptr);
    char stime[20];
    std::strftime(stime, sizeof(stime), "%Y-%m-%d %H:%M:%S", std::localtime(&time));
    return std::string(stime);
}

#define LOG std::cout << getTimestamp() << " : "
#define LOG_ERR std::cerr << getTimestamp() << " : "

class Function4D {
   public:
    virtual double operator()(double t, double x, double y, double z) const = 0;
};

class Function3D {
   public:
    virtual double operator()(double x, double y, double z) const = 0;
};

class Grid3D {
    /* Raw data */
    double *raw;

    const int _N;

   public:
    /* Grid size is (N + 1) ^ 3 */
    const int N;
    const int size;

    Grid3D(int N) : N(N), _N(N + 1), size(_N * _N * _N) {
        raw = new double[size];
        LOG << "Created raw array of size " << size << endl;
    }

    ~Grid3D() {
        if (raw) {
            delete raw;
            LOG << "Deleted raw array of size " << size << endl;
        }
    }

    double &operator()(int i, int j, int k) {
        int idx = i * (_N * _N) + j * _N + k;
        return raw[idx];
    }

    double operator()(int i, int j, int k) const {
        int idx = i * (_N * _N) + j * _N + k;
        return raw[idx];
    }

    void writeToFile(std::ofstream &outFile) const {
        outFile.write((char *)raw, size * sizeof(double));
    }
};

class Solver {
    const Function4D *const u;
    const Function3D *const phi;

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

   public:
    Solver(double T, double L_x, double L_y, double L_z, int N, int K, Function4D *u, Function3D *phi)
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

    void init_0(Grid3D &grid) {
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

    void init_1(Grid3D &grid) {
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

    double laplasian(const Grid3D &g, int i, int j, int k) {
        double center = g(i, j, k);
        return (g(i - 1, j, k) - 2.0 * center + g(i + 1, j, k)) / h_x +
               (g(i, j - 1, k) - 2.0 * center + g(i, j + 1, k)) / h_y +
               (g(i, j, k - 1) - 2.0 * center + g(i, j, k + 1)) / h_z;
    }

    /** makeStep fills n-th layer of grid. It depends on two previous layers. */
    void makeStep(const int n, Grid3D &layer, const Grid3D &previous_1, const Grid3D &previous_2) {
        if (n < 2 || n >= K) {
            LOG_ERR << "Parameter n in makeStep must be between 2 and T. Actual value: " << n << endl;
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
            LOG << "Step " << n + 1 << " completed" << endl;
        }
    }

    Grid3D *solve() {
        LOG << "Start solve()\n";
        Grid3D *grids[3] = {new Grid3D(N), new Grid3D(N), new Grid3D(N)};
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

    void fillByF(Grid3D &grid, Function4D *f, int n) {
        for (int i = 0; i <= grid.N; ++i) {
            for (int j = 0; j <= grid.N; ++j) {
                for (int k = 0; k <= grid.N; ++k) {
                    grid(i, j, k) = (*f)(
                        tau * n,
                        h_x * i,
                        h_y * j,
                        h_z * k);
                }
            }
        }
    }
};

class U : public Function4D {
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

class Phi : public Function3D {
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
    if (argc <= 8) {
        LOG_ERR << "Usage: prog gt|num L_x L_y L_z T N K out_file" << endl;
        return 0;
    }

    std::ofstream outFile(argv[8], std::ios::out | std::ios::binary);
    LOG << "Output file created\n";

    double L_x = atof(argv[2]);
    double L_y = atof(argv[3]);
    double L_z = atof(argv[4]);
    double T = atof(argv[5]);
    int N = atoi(argv[6]);
    int K = atoi(argv[7]);
    LOG << "Papameters parsed succesfully\n";

    Phi phi(L_x, L_y, L_z);
    U u(L_x, L_y, L_z);
    LOG << "Phi and U created\n";

    Solver solver(K, L_x, L_y, L_z, N, K, &u, &phi);
    LOG << "Solver created\n";
    LOG << "Initialization successfully completed\n";
    Grid3D *grid;

    if (strcmp(argv[1], "num") == 0) {
        grid = solver.solve();
        LOG << "Solving complete\n";
    } else {
        grid = new Grid3D(N);
        solver.fillByF(*grid, &u, K);
    }

    LOG << "Writing result to file\n";
    grid->writeToFile(outFile);
    outFile.close();
    LOG << "Result written\n";

    delete grid;
    return 0;
}