#pragma once

#include "functions.h"
#include "Grid3D.h"

class MathSolver {
public:
    MathSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi);

    void init_0(Grid3D &grid, int start_i, int start_j, int start_k) const;

    void init_1(Grid3D &grid, int start_i, int start_j, int start_k) const;

    /** Fills n-th grid of grid. It depends on two previous layers. */
    void fillInnerNodes(Grid3D &grid, const Grid3D &previous_1, const Grid3D &previous_2) const;

    double laplacian(const Grid3D &g, int i, int j, int k) const;

    void fillByU(Grid3D &grid, int n) const;

    void C_norm(Grid3D &one, Grid3D &another) const;

private:
    const U u;
    const Phi phi;

    const int N;
    const int K;

    /* Grid steps */
    const double h_x;
    const double h_y;
    const double h_z;
    const double tau;
};