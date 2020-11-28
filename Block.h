#pragma once

#include <vector>
#include <mpi.h>
#include "Grid3D.h"
#include "MathSolver.h"

const int N_GRIDS = 3;

class Block {
public:
    static int calcBlockSize(int N, int I, int splits);

    Block(MathSolver *solver, int splits_X, int splits_Y, int splits_Z, int N);

    void makeStep();

    const Grid3D *getCurrentState() const;

private:
    int block_coords[3];
    int block_size[3];
    int start[3];
    int n_splits[3];   // Number of blocks
    int blockId;
    bool isPeriodicalCondition[3];
    const MathSolver *solver;

    /* Mutable state */
    int iteration;
    std::vector<Grid3D> grids;

    void syncWithNeighbors();

    int getNeighborId(int axis, int direction) const;
};
