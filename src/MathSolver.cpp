#include "MathSolver.h"

MathSolver::MathSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi)
        : _N(N),
          K(K),
          u(u),
          phi(phi),
          h_x(L_x / N),
          h_y(L_y / N),
          h_z(L_z / N),
          tau(T / K) {}


std::ostream &operator<<(std::ostream &out, const MathSolver *solver) {
    return out << "MathSolver: "
               << "N = " << solver->_N << ", "
               << "K = " << solver->K << ", "
               << "h_x = " << solver->h_x << ", "
               << "h_y = " << solver->h_y << ", "
               << "h_z = " << solver->h_z << ", "
               << "tau = " << solver->tau;
}
