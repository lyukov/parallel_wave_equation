#pragma once

#include "Grid3D.h"
#include "functions.h"

class MathSolver {
public:
    MathSolver(double T, double L_x, double L_y, double L_z, int N, int K, U u, Phi phi);

    virtual ~MathSolver() {}

    virtual void init_0(int start_i, int start_j, int start_k) = 0;

    virtual void init_1() = 0;

    virtual void makeStepForInnerNodes(int n) = 0;

    virtual void updateGroundTruth(int n, int start_i, int start_j, int start_k) = 0;

    virtual double maxAbsoluteErrorInner(int n) = 0;

    virtual double sumSquaredErrorInner(int n) = 0;

    virtual std::vector<double> getSlice(int n, int index, int axis) = 0;

    virtual int getSliceSize(int axis) = 0;

    virtual void setSlice(int n, int index, int axis, std::vector<double> &slice) = 0;

    virtual void setZeros(int n, int index, int axis) = 0;

    virtual double maxGroundTruth();

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
