#pragma once
#include <fstream>

class Grid3D {
    double *raw;
    const int _N1, _N2, _N3;

   public:
    const int size;

    Grid3D(int N);

    Grid3D(int N1, int N2, int N3);

    ~Grid3D();

    double &operator()(int i, int j, int k);

    double operator()(int i, int j, int k) const;

    void writeToFile(std::ofstream &outFile) const;
};