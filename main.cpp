#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <string>

#include "Grid3D.h"
#include "Solver.h"
#include "functions.h"
#include "log.h"

using std::endl;

std::string getTimestamp() {
    std::time_t time = std::time(NULL);
    char stime[20];
    std::strftime(stime, sizeof(stime), "%Y-%m-%d %H:%M:%S", std::localtime(&time));
    return std::string(stime);
}

int main(int argc, char **argv) {
    if (argc <= 8) {
        LOG_ERR << "Usage: prog gt|num L_x L_y L_z T N K out_file" << endl;
        return 0;
    }

    char *target = argv[1];
    double L_x = atof(argv[2]);
    double L_y = atof(argv[3]);
    double L_z = atof(argv[4]);
    double T = atof(argv[5]);
    int N = atoi(argv[6]);
    int K = atoi(argv[7]);
    char *outFileName = argv[8];
    LOG << "Papameters parsed succesfully\n";

    std::ofstream outFile(outFileName, std::ios::out | std::ios::binary);
    LOG << "Output file created\n";

    Phi phi(L_x, L_y, L_z);
    U u(L_x, L_y, L_z);
    LOG << "Phi and U created\n";

    Solver solver(T, L_x, L_y, L_z, N, K, u, phi);
    LOG << "Solver created\n";
    LOG << "Initialization successfully completed\n";
    Grid3D *grid;

    if (strcmp(target, "num") == 0) {
        grid = solver.solve();
        LOG << "Solving complete\n";
    } else if (strcmp(target, "gt") == 0) {
        grid = new Grid3D(N);
        solver.fillByU(*grid, K);
    } else {
        LOG_ERR << "Unknown target" << endl;
        return 1;
    }

    LOG << "Writing result to file\n";
    grid->writeToFile(outFile);
    outFile.close();
    LOG << "Result written\n";

    delete grid;
    return 0;
}