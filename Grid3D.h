#pragma once

#include <fstream>

/**
 * Three-dimensional tensor, that stores values [0, ..., N]
 */
class Grid3D {
    double *raw;
    const int _N1;
    const int _N2;
    const int _N3;
    const int offset_1;
    const int offset_2;
    const int offset_3;

    int _index(int i, int j, int k) const;

    void init();

public:
    const int size;

    explicit Grid3D(int N);

    Grid3D(int N1, int N2, int N3);

    Grid3D(int N1, int N2, int N3, int offset_1, int offset_2, int offset_3);

    ~Grid3D();

    double &operator()(int i, int j, int k);

    double operator()(int i, int j, int k) const;

    void writeToFile(std::ofstream &outFile) const;
};