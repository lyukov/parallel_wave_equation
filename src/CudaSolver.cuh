#pragma once

#include "Grid3D.h"

void makeStepWithCuda(Grid3D &grid, const Grid3D &previous_1, const Grid3D &previous_2,
                      double h_x, double h_y, double h_z, double sqr_tau);