#pragma once

#include "Grid3D.h"
#include "functions.h"
#include "MathSolver.h"

class CudaSolver : public MathSolver {
public:
    CudaSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi);

    void makeStepForInnerNodes(Grid3D &grid, Grid3D &previous_1, Grid3D &previous_2);

    void fillByGroundTruth(Grid3D &grid, int n, int start_i, int start_j, int start_k);

    double maxAbsoluteErrorInner(const Grid3D &grid, const Grid3D &another);
};