#pragma once

#include "functions.h"
#include "Grid3D.h"
#include "utils.h"
#include <iostream>

class MathSolver {
public:
    MathSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi);

    void init_0(Grid3D &grid, int start_i, int start_j, int start_k) const;

    void init_1(Grid3D &grid, Grid3D &previous) const;

    /** Fills n-th grid of grid. It depends on two previous layers. */
    void makeStepForInnerNodes(Grid3D &grid, const Grid3D &previous_1, const Grid3D &previous_2) const;

    void fillByGroundTruth(Grid3D &grid, int n, int start_i, int start_j, int start_k) const;

    double maxAbsoluteErrorInner(const Grid3D &grid, const Grid3D &another) const;

    double sumSquaredErrorInner(const Grid3D &grid, const Grid3D &another) const;

    double maxRelativeErrorInner(const Grid3D &grid, const Grid3D &another) const;

private:
    const U u;
    const Phi phi;

    const int _N;
    const int K;

    /* Grid steps */
    const double h_x;
    const double h_y;
    const double h_z;
    const double tau;

    friend std::ostream &operator<<(std::ostream &out, const MathSolver &solver);

    double laplacian(const Grid3D &g, int i, int j, int k) const;

    double laplacian(const Phi &phi, double x, double y, double z, double h) const;
};

std::ostream &operator<<(std::ostream &out, const MathSolver &solver);