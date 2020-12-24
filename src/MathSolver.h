#pragma once

#include "functions.h"
#include "Grid3D.h"
#include "utils.h"
#include <iostream>

class MathSolver {
public:
    MathSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi);

    virtual void init_0(Grid3D &grid, int start_i, int start_j, int start_k);

    virtual void init_1(Grid3D &grid, Grid3D &previous);

    virtual void makeStepForInnerNodes(Grid3D &grid, Grid3D &previous_1, Grid3D &previous_2);

    virtual void fillByGroundTruth(Grid3D &grid, int n, int start_i, int start_j, int start_k);

    virtual double maxAbsoluteErrorInner(const Grid3D &grid, const Grid3D &another);

    virtual double sumSquaredErrorInner(const Grid3D &grid, const Grid3D &another);

    virtual double maxRelativeErrorInner(const Grid3D &grid, const Grid3D &another);

protected:
    const U u;
    const Phi phi;

    const int _N;
    const int K;

    /* Grid steps */
    const double h_x;
    const double h_y;
    const double h_z;
    const double tau;

    friend std::ostream &operator<<(std::ostream &out, const MathSolver *solver);

    double laplacian(const Grid3D &g, int i, int j, int k) const;
};

std::ostream &operator<<(std::ostream &out, const MathSolver *solver);