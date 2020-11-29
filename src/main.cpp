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
    if (argc <= 9) {
        LOG_ERR << "argc = " << argc << endl;
        LOG_ERR << "Usage: prog L_x L_y L_z T N K splits_X splits_Y splits_Z" << endl;
        return 0;
    }

    double L_x = atof(argv[1]);
    double L_y = atof(argv[2]);
    double L_z = atof(argv[3]);
    double T = atof(argv[4]);
    int N = atoi(argv[5]);
    int K = atoi(argv[6]);
    int splits_X = atoi(argv[7]);
    int splits_Y = atoi(argv[8]);
    int splits_Z = atoi(argv[9]);
    LOG << "Papameters parsed succesfully\n";

    MPIProxy mpiProxy(&argc, &argv);
    int nProcessors = mpiProxy.getNumOfProcessors();
    assert(nProcessors == splits_X * splits_Y * splits_Z);
    LOG << "MPI Proxy created. Rank: " << mpiProxy.getRank() << ". Processors: " << nProcessors << endl;

//    std::ofstream outFile(outFileName, std::ios::out | std::ios::binary);
//    LOG << "Output file created\n";

    Phi phi(L_x, L_y, L_z);
    U u(L_x, L_y, L_z);
    LOG << "Phi and U created" << endl;

    MathSolver solver(T, L_x, L_y, L_z, N, K, u, phi);
    LOG << "MathSolver created" << endl;
    Block block(&mpiProxy, &solver, splits_X, splits_Y, splits_Z, N);
    LOG << "Block created" << endl;
    LOG << "Initialization successfully completed" << endl;

    Grid3D groundTruth(block.shape[0], block.shape[1], block.shape[2]);

    for (int iter = 0; iter < K; ++iter) {
        block.makeStep();
        solver.fillByU(
                groundTruth,
                iter,
                block.start[0] - 1,
                block.start[1] - 1,
                block.start[2] - 1
        );
        block.printError(groundTruth);
    }

    return 0;
}