#pragma once

#include "Grid3D.h"

void makeStepWithCuda(Grid3D &grid, Grid3D &previous_1, Grid3D &previous_2,
                      double h_x, double h_y, double h_z, double sqr_tau);