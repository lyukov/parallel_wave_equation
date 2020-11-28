#include "Block.h"
#include "utils.h"

int Block::calcBlockSize(int N, int I, int splits) {
    const int size = N / splits;
    return (I == splits - 1) ?         // if last block
           N - size * (splits - 1) :   // get all remaining size
           size;
}

Block::Block(MathSolver *solver, int splits_X, int splits_Y, int splits_Z, int N) : solver(solver) {
    MPI_Comm_rank(MPI_COMM_WORLD, &blockId);

    iteration = 0;
    block_coords[0] = blockId / (splits_Y * splits_Z);
    block_coords[1] = (blockId - splits_Y * splits_Z * block_coords[0]) / splits_Z;
    block_coords[2] = blockId % splits_Z;
    block_size[0] = calcBlockSize(N, block_coords[0], splits_X);
    block_size[1] = calcBlockSize(N, block_coords[1], splits_Y);
    block_size[2] = calcBlockSize(N, block_coords[2], splits_Z);
    n_splits[0] = splits_X;
    n_splits[1] = splits_Y;
    n_splits[2] = splits_Z;
    start[0] = block_coords[0] * block_size[0];
    start[1] = block_coords[1] * block_size[1];
    start[2] = block_coords[2] * block_size[2];
    isPeriodicalCondition[0] = true;
    isPeriodicalCondition[1] = false;
    isPeriodicalCondition[2] = true;

    for (int i = 0; i < 3; ++i) {
        grids.push_back(
                Grid3D(block_size[0] + 2, block_size[1] + 2, block_size[2] + 2) // Add padding
        );
    }
}

const Grid3D &Block::getCurrentState() const {
    return grids[(iteration + N_GRIDS - 1) % N_GRIDS];
}

void Block::makeStep() {
    if (iteration == 0) {
        solver->init_0(grids[iteration % N_GRIDS], start[0], start[1], start[2]);
    } else if (iteration == 1) {
        solver->init_1(grids[iteration % N_GRIDS], start[0], start[1], start[2]);
    } else {
        solver->fillInnerNodes(
                grids[iteration % N_GRIDS],
                grids[(iteration - 1) % N_GRIDS],
                grids[(iteration - 2) % N_GRIDS]
        );
        syncWithNeighbors();
    }
    iteration++;
    if ((iteration % 50) == 0) {
        LOG << "Iteration " << iteration << " completed" << endl;
    }
}

void Block::syncWithNeighbors() {
    MPI_Request dummyRequest;
    Grid3D &grid = grids[iteration % N_GRIDS];

    for (int axis = 0; axis < 3; ++axis) {
        for (int direction = -1; direction <= 1; direction += 2) {
            int neighborId = getNeighborId(axis, direction);
            if (neighborId == -1) continue;
            int index = direction == -1 ? 0 : block_size[axis] + 1;
            Slice slice = grid.getSlice(index, axis);
            MPI_Isend(slice.data(), slice.size(), MPI_DOUBLE,
                      neighborId, 0, MPI_COMM_WORLD, &dummyRequest);
        }
    }

    MPI_Status dummyStatus;
    for (int axis = 0; axis < 3; ++axis) {
        for (int direction = -1; direction <= 1; direction += 2) {
            int neighborId = getNeighborId(axis, direction);
            if (neighborId == -1) continue;
            int index = direction == -1 ? 0 : block_size[axis] + 1;
            Slice slice(grid.getSliceSize(axis));
            MPI_Recv(slice.data(), slice.size(), MPI_DOUBLE,
                     neighborId, 0, MPI_COMM_WORLD, &dummyStatus);
            grid.setSlice(index, axis, slice);
        }
    }
}

int Block::getNeighborId(int axis, int direction) const {
    int neighbor[3] = {block_coords[0], block_coords[1], block_coords[2]};
    neighbor[axis] = neighbor[axis] + direction;
    if (neighbor[axis] < 0 || neighbor[axis] >= n_splits[axis]) {
        if (isPeriodicalCondition[axis]) return -1;
        neighbor[axis] = (neighbor[axis] + n_splits[axis]) % n_splits[axis];
    }
    return neighbor[0] * n_splits[1] * n_splits[2] +
           neighbor[1] * n_splits[2] +
           neighbor[2];
}

void Block::printError(Grid3D &groundTruth) const {
    double error = solver->C_norm_inner(getCurrentState(), groundTruth);
    double wholeError;
    MPI_Reduce(&error, &wholeError, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (blockId == 0) {
        LOG << "Iteration: " << iteration << ". Error: " << wholeError << endl;
    }
}
