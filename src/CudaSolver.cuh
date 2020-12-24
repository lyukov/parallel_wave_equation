#pragma once

#include "Grid3D.h"
#include "functions.h"

void makeStepWithCuda(Grid3D &grid, Grid3D &previous_1, Grid3D &previous_2,
                      double h_x, double h_y, double h_z, double sqr_tau);

void fillByGtWithCuda(Grid3D &grid, U u, int n, double tau, int start_i, int start_j, int start_k);