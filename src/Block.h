#pragma once

#include <vector>
#include "MPIProxy.h"
#include "Grid3D.h"
#include "CpuSolver.h"

class Block {
public:
    int shape[3];
    int start[3];

    static int calcBlockSize(int N, int I, int splits);

    Block(MPIProxy *mpi, int splits_X, int splits_Y, int splits_Z, int N);

    void makeStep();

    double printError();

    void setSolver(MathSolver* solver);

private:
    int N;
    int block_coords[3];
    int n_splits[3];   // Number of blocks
    bool isPeriodicalCondition[3];
    MathSolver *solver;
    const MPIProxy *mpi;

    /* Mutable state */
    int iteration;

    void syncWithNeighbors();

    int getNeighborId(int axis, int direction) const;

    void sendToNeighbor(int axis, int direction, std::vector<double> &buf);

    void receiveFromNeighbor(int axis, int direction);
};
