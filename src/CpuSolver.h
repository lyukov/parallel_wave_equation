#pragma once

#include "functions.h"
#include "Grid3D.h"
#include "utils.h"
#include "MathSolver.h"
#include <iostream>

class CpuSolver : public MathSolver {
public:
    CpuSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi);

    void init_0(Grid3D &grid, int start_i, int start_j, int start_k) override;

    void init_1(Grid3D &grid, Grid3D &previous) override;

    void makeStepForInnerNodes(Grid3D &grid, Grid3D &previous_1, Grid3D &previous_2) override;

    void fillByGroundTruth(Grid3D &grid, int n, int start_i, int start_j, int start_k) override;

    double maxAbsoluteErrorInner(Grid3D &grid, Grid3D &another) override;

    double sumSquaredErrorInner(Grid3D &grid, Grid3D &another) override;

protected:
    double laplacian(const Grid3D &g, int i, int j, int k) const;
};
