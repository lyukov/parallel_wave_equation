#pragma once

#include "functions.h"
#include "Grid3D.h"
#include "utils.h"
#include "MathSolver.h"
#include <iostream>

class CpuSolver : public MathSolver {
public:
    CpuSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi,
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

    double maxGroundTruth() override;

protected:
    const int N_GRIDS = 3;

    double laplacian(const Grid3D &g, int i, int j, int k) const;

    std::vector<Grid3D> grids;

    Grid3D groundTruth;

    Grid3D &getCurrentState(int n);
};
