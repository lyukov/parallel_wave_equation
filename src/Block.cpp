#include "Block.h"
#include "utils.h"

int Block::calcBlockSize(int N, int I, int splits) {
    int size = N / splits;
    size = (I == splits - 1) ?         // if last block
           N - size * (splits - 1) :   // get all remaining size
           size;
    return size + 2;                   // padding
}

Block::Block(
        MPIProxy *mpi,
        int splits_X, int splits_Y, int splits_Z,
        int N
) : solver(0), mpi(mpi), N(N) {
    int blockId = mpi->getRank();
    iteration = 0;
    isPeriodicalCondition[0] = true;
    isPeriodicalCondition[1] = false;
    isPeriodicalCondition[2] = true;
    block_coords[0] = blockId / (splits_Y * splits_Z);
    block_coords[1] = (blockId - splits_Y * splits_Z * block_coords[0]) / splits_Z;
    block_coords[2] = blockId % splits_Z;
    shape[0] = calcBlockSize(N, block_coords[0], splits_X);
    shape[1] = calcBlockSize(N, block_coords[1], splits_Y);
    shape[2] = calcBlockSize(N, block_coords[2], splits_Z);
    n_splits[0] = splits_X;
    n_splits[1] = splits_Y;
    n_splits[2] = splits_Z;
    start[0] = block_coords[0] * (N / n_splits[0]);
    start[1] = block_coords[1] * (N / n_splits[1]);
    start[2] = block_coords[2] * (N / n_splits[2]);
}

void Block::makeStep() {
    if (iteration == 0) {
        solver->init_0(start[0] - 1, start[1] - 1, start[2] - 1);
    } else if (iteration == 1) {
        solver->init_1();
    } else {
        solver->makeStepForInnerNodes(iteration);
    }
    syncWithNeighbors();
    iteration++;
}

void Block::syncWithNeighbors() {
    std::vector<double> buf;
    for (int axis = 0; axis < 3; ++axis) {
        for (int direction = -1; direction <= 1; direction += 2) {
            sendToNeighbor(axis, direction, buf);
            receiveFromNeighbor(axis, -direction);
            mpi->barrier();
        }
    }
}

void Block::sendToNeighbor(int axis, int direction, std::vector<double> &buf) {
    int neighborId = getNeighborId(axis, direction);
    if (neighborId != -1) {
        int index = direction == -1 ? 1 : shape[axis] - 2;
        buf = solver->getSlice(iteration, index, axis);
        mpi->sendVector(buf, neighborId);
    }
}

void Block::receiveFromNeighbor(int axis, int direction) {
    int neighborId = getNeighborId(axis, direction);
    if (neighborId != -1) {
        int sliceSize = solver->getSliceSize(axis);
        std::vector<double> slice = mpi->receiveVector(sliceSize, neighborId);
        int index = direction == -1 ? 0 : shape[axis] - 1;
        solver->setSlice(iteration, index, axis, slice);
    } else {
        // Zero border conditions
        int index = direction == -1 ? 1 : shape[axis] - 1;
        solver->setZeros(iteration, index, axis);
    }
}

int Block::getNeighborId(int axis, int direction) const {
    int neighbor[3] = {block_coords[0], block_coords[1], block_coords[2]};
    neighbor[axis] = neighbor[axis] + direction;
    if (neighbor[axis] < 0 || neighbor[axis] >= n_splits[axis]) {
        if (!isPeriodicalCondition[axis]) {
            return -1;
        }
        neighbor[axis] = (neighbor[axis] + n_splits[axis]) % n_splits[axis];
    }
    return neighbor[0] * n_splits[1] * n_splits[2] +
           neighbor[1] * n_splits[2] +
           neighbor[2];
}

double Block::printError(Grid3D &groundTruth) {
    //if (iteration % 10) return;
    double absoluteError = mpi->maxOverAll(
            solver->maxAbsoluteErrorInner(iteration - 1)
    );
    mpi->barrier();
    double sumSquaredError = mpi->sumOverAll(
            solver->sumSquaredErrorInner(iteration - 1)
    );
    double mse = sumSquaredError / (N * N * N);
    mpi->barrier();
    double maxGt = mpi->maxOverAll(
            max(groundTruth.getFlatten())
    );
    if (mpi->isMainProcess()) {
        LOG << "Iteration: " << iteration
            << ". \tC1: " << absoluteError
            << " \tMSE: " << mse
            << " \tMax GT: " << maxGt << endl;
    }
    return absoluteError;
}

void Block::setSolver(MathSolver *solver) {
    this->solver = solver;
}
