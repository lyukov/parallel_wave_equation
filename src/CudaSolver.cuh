#pragma once

#include "Grid3D.h"
#include "functions.h"
#include "MathSolver.h"

class CudaSolver : public MathSolver {
public:
    CudaSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi,
               int shapeX, int shapeY, int shapeZ);

    void init_0(int start_i, int start_j, int start_k) override;

    void init_1() override;

    void makeStepForInnerNodes(int n) override;

    void updateGroundTruth(int n, int start_i, int start_j, int start_k) override;

    double maxAbsoluteErrorInner(int n) override;

    double sumSquaredErrorInner(int n) override;

    std::vector<double> getSlice(int n, int index, int axis) override;

    int getSliceSize(int axis) override;

    void setSlice(int n, int index, int axis, std::vector<double> &slice) override;

    void setZeros(int n, int index, int axis) override;

private:
    const int N_GRIDS = 3;

    std::vector<double *> d_grids;
    double *d_groundTruth;
    double *d_errorC1, *d_errorMSE;
    Grid3D grid3D;

    const unsigned long sizeInBytes;
    const unsigned long flatSize;
    const int blockSizeFull;
    const int blockSizeInner;
    const int gridSizeFull;
    const int gridSizeInner;
};