#pragma once

#include "functions.h"
#include "Grid3D.h"

class Solver {
    const Function4D *const u;
    const Function3D *const phi;

    const int N;
    const int K;

    /* Parameters of continuous problem */
    const double T;
    const double L_x;
    const double L_y;
    const double L_z;

    /* Grid steps */
    const double h_x;
    const double h_y;
    const double h_z;
    const double tau;

public:
    Solver(double T, double L_x, double L_y, double L_z, int N, int K, Function4D *u, Function3D *phi);

    void init_0(Grid3D &grid);

    void init_1(Grid3D &grid);

    double laplasian(const Grid3D &g, int i, int j, int k) const;

    /** Fills n-th layer of grid. It depends on two previous layers. */
    void makeStep(const int n, Grid3D &layer, const Grid3D &previous_1, const Grid3D &previous_2) const;

    Grid3D *solve();

    void fillByF(Grid3D &grid, Function4D *f, int n);
};