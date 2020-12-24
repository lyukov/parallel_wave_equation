#pragma once

#include "Grid3D.h"
#include "functions.h"
#include "CpuSolver.h"

class CudaSolver : public CpuSolver {
public:
    CudaSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi, Grid3D &grid);

    void init_0(Grid3D &grid, int start_i, int start_j, int start_k);

    void init_1(Grid3D &grid, Grid3D &previous);

    void makeStepForInnerNodes(Grid3D &grid, Grid3D &previous_1, Grid3D &previous_2);

    void fillByGroundTruth(Grid3D &grid, int n, int start_i, int start_j, int start_k);

    double maxAbsoluteErrorInner(Grid3D &grid, Grid3D &another);

    double sumSquaredErrorInner(Grid3D &grid, Grid3D &another);

private:
    const unsigned long sizeInBytes;
    const unsigned long flatSize;
    const int blockSizeFull;
    const int blockSizeInner;
    const int gridSizeFull;
    const int gridSizeInner;
};