#include <ctime>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#include "CudaSolver.h"
#include "CpuSolver.h"
#include "functions.h"
#include "utils.h"
#include "MPIProxy.h"
#include "Block.h"

std::string getTimestamp() {
    std::time_t time = std::time(NULL);
    char stime[20];
    std::strftime(stime, sizeof(stime), "%Y-%m-%d %H:%M:%S", std::localtime(&time));
    return std::string(stime);
}

int main(int argc, char **argv) {
    if (argc <= 8) {
        LOG_ERR << "argc = " << argc << endl;
        LOG_ERR << "Usage: wave L T N K splits_X splits_Y splits_Z label" << endl;
        return 0;
    }

    double L = atof(argv[1]);
    double T = atof(argv[2]);
    int N = atoi(argv[3]);
    int K = atoi(argv[4]);
    int splits_X = atoi(argv[5]);
    int splits_Y = atoi(argv[6]);
    int splits_Z = atoi(argv[7]);
    std::string label(argv[8]);
    double L_x = L, L_y = L, L_z = L;

    std::stringstream csvOut;

    MPIProxy mpi(&argc, &argv);
    double startTime = mpi.time();
    int nProcessors = mpi.getNumOfProcessors();
    if (nProcessors != splits_X * splits_Y * splits_Z) {
        LOG_ERR << "Incorrect num of processors" << endl;
        return 1;
    }
    // LOG_DEBUG << "MPI Proxy created. Rank: " << mpi.getRank() << ". Processors: " << nProcessors << endl;

    Phi phi(L_x, L_y, L_z);
    U u(L_x, L_y, L_z);

    Block block(&mpi, splits_X, splits_Y, splits_Z, N);

    Grid3D groundTruth(block.shape[0], block.shape[1], block.shape[2]);

    MathSolver *solver = new CudaSolver(
            T, L_x, L_y, L_z, N, K, u, phi,
            block.shape[0], block.shape[1], block.shape[2]
    );
    if (mpi.isMainProcess()) {
        csvOut << label << TAB << nProcessors
               << TAB << L << TAB << T << TAB << N << TAB << K
               << TAB << L / N << TAB << T / K;
        LOG << "Num of processors: " << nProcessors << endl;
        LOG << solver << endl;
    }
    block.setSolver(solver);

    double error;
    for (int iter = 0; iter <= K; ++iter) {
        block.makeStep();
        solver->updateGroundTruth(
                iter,
                block.start[0] - 1,
                block.start[1] - 1,
                block.start[2] - 1
        );
        error = block.printError();
    }

//    /* Save result as binary file */
//    if (label == "dump") {
//        std::stringstream ss;
//        ss << "block_" << mpi.getRank();
//        block.getCurrentState().writeToFile(ss.str());
//        ss << "_gt";
//        groundTruth.writeToFile(ss.str());
//    }

    mpi.barrier();
    double duration = mpi.time() - startTime;
    if (mpi.isMainProcess()) {
        csvOut << TAB << error << TAB << duration;
        std::ofstream csvStats("stats.csv", std::ios_base::app);
        csvStats << csvOut.str() << endl;
        csvStats.close();
        LOG << "All processes finished. Elapsed time (s): " << duration << endl;
    }

    return 0;
}