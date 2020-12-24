#pragma once

#include "Grid3D.h"
#include "functions.h"

class MathSolver {
public:
    MathSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi);

    virtual ~MathSolver() {}

    virtual void init_0(Grid3D &grid, int start_i, int start_j, int start_k) = 0;

    virtual void init_1(Grid3D &grid, Grid3D &previous) = 0;

    virtual void makeStepForInnerNodes(Grid3D &grid, Grid3D &previous_1, Grid3D &previous_2) = 0;

    virtual void fillByGroundTruth(Grid3D &grid, int n, int start_i, int start_j, int start_k) = 0;

    virtual double maxAbsoluteErrorInner(Grid3D &grid, Grid3D &another) = 0;

    virtual double sumSquaredErrorInner(Grid3D &grid, Grid3D &another) = 0;

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

};

std::ostream &operator<<(std::ostream &out, const MathSolver *solver);
