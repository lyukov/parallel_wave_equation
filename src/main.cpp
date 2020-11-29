#include <ctime>
#include <iostream>
#include <string>

#include "MathSolver.h"
#include "functions.h"
#include "utils.h"
#include "MPIProxy.h"
#include "Block.h"

using std::endl;

std::string getTimestamp() {
    std::time_t time = std::time(NULL);
    char stime[20];
    std::strftime(stime, sizeof(stime), "%Y-%m-%d %H:%M:%S", std::localtime(&time));
    return std::string(stime);
}

int main(int argc, char **argv) {
    if (argc <= 7) {
        LOG_ERR << "argc = " << argc << endl;
        LOG_ERR << "Usage: prog L T N K splits_X splits_Y splits_Z" << endl;
        return 0;
    }

    double L = atof(argv[1]);
    double T = atof(argv[2]);
    int N = atoi(argv[3]);
    int K = atoi(argv[4]);
    int splits_X = atoi(argv[5]);
    int splits_Y = atoi(argv[6]);
    int splits_Z = atoi(argv[7]);
    double L_x = L, L_y = L, L_z = L;
    LOG_DEBUG << "Papameters parsed succesfully\n";

    MPIProxy mpi(&argc, &argv);
    double startTime = mpi.time();
    int nProcessors = mpi.getNumOfProcessors();
    assert(nProcessors == splits_X * splits_Y * splits_Z);
    LOG_DEBUG << "MPI Proxy created. Rank: " << mpi.getRank() << ". Processors: " << nProcessors << endl;

//    std::ofstream outFile(outFileName, std::ios::out | std::ios::binary);
//    LOG << "Output file created\n";

    Phi phi(L_x, L_y, L_z);
    U u(L_x, L_y, L_z);

    MathSolver solver(T, L_x, L_y, L_z, N, K, u, phi);
    if (mpi.isMainProcess()) {
        LOG << solver << endl;
    }
    Block block(&mpi, &solver, splits_X, splits_Y, splits_Z, N);
    LOG_DEBUG << "Block created" << endl;
    LOG_DEBUG << "Initialization successfully completed" << endl;

    Grid3D groundTruth(block.shape[0], block.shape[1], block.shape[2]);

    for (int iter = 0; iter <= K; ++iter) {
        block.makeStep(iter != K);
        solver.fillByGroundTruth(
                groundTruth,
                iter,
                block.start[0] - 1,
                block.start[1] - 1,
                block.start[2] - 1
        );
        block.printError(groundTruth);
    }

    mpi.barrier();
    double duration = mpi.time() - startTime;
    if (mpi.isMainProcess()) {
        LOG << "All processes finished. Elapsed time (s): " << duration << endl;
    }

    return 0;
}