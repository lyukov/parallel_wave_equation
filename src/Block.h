#pragma once

#include <vector>
#include "MPIProxy.h"
#include "Grid3D.h"
#include "MathSolver.h"

const int N_GRIDS = 3;

class Block {
public:
    int shape[3];
    int start[3];

    static int calcBlockSize(int N, int I, int splits);

    Block(MPIProxy *mpi, MathSolver *solver, int splits_X, int splits_Y, int splits_Z, int N);

    void makeStep(bool shareBorders = true);

    const Grid3D &getCurrentState() const;

    void printError(Grid3D &groundTruth) const;

private:
    int block_coords[3];
    int n_splits[3];   // Number of blocks
    bool isPeriodicalCondition[3];
    const MathSolver *solver;
    const MPIProxy *mpi;

    /* Mutable state */
    int iteration;
    std::vector<Grid3D> grids;

    void syncWithNeighbors();

    int getNeighborId(int axis, int direction) const;

    void sendToNeighbors(int axis, int direction);

    void receiveFromNeighbors(int axis, int direction);
};
