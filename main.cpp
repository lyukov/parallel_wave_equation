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

U::U(double L_x, double L_y, double L_z)
        : L_x(L_x), L_y(L_y), L_z(L_z),
          a_t(M_PI * sqrt(4.0 / (L_x * L_x) + 1.0 / (L_y * L_y) + 1.0 / (L_z * L_z))) {}

double U::operator()(double t, double x, double y, double z) const {
    return sin(M_2_PI * x / L_x) * sin(M_PI * y / L_y) * sin(M_2_PI * z / L_z) * cos(a_t * t);
}

Phi::Phi(double L_x, double L_y, double L_z) : L_x(L_x), L_y(L_y), L_z(L_z) {}

double Phi::operator()(double x, double y, double z) const {
    return sin(M_2_PI * x / L_x) * sin(M_PI * y / L_y) * sin(M_2_PI * z / L_z);
}

int main(int argc, char **argv) {
    if (argc <= 8) {
        LOG_ERR << "Usage: prog gt|num L_x L_y L_z T N K out_file" << endl;
        return 0;
    }

    std::ofstream outFile(argv[8], std::ios::out | std::ios::binary);
    LOG << "Output file created\n";

    double L_x = atof(argv[2]);
    double L_y = atof(argv[3]);
    double L_z = atof(argv[4]);
    double T = atof(argv[5]);
    int N = atoi(argv[6]);
    int K = atoi(argv[7]);
    LOG << "Papameters parsed succesfully\n";

    Phi phi(L_x, L_y, L_z);
    U u(L_x, L_y, L_z);
    LOG << "Phi and U created\n";

    Solver solver(T, L_x, L_y, L_z, N, K, &u, &phi);
    LOG << "Solver created\n";
    LOG << "Initialization successfully completed\n";
    Grid3D *grid;

    if (strcmp(argv[1], "num") == 0) {
        grid = solver.solve();
        LOG << "Solving complete\n";
    } else {
        grid = new Grid3D(N);
        solver.fillByF(*grid, &u, K);
    }

    LOG << "Writing result to file\n";
    grid->writeToFile(outFile);
    outFile.close();
    LOG << "Result written\n";

    delete grid;
    return 0;
}