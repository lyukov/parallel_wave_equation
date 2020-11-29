#pragma once

#include <fstream>
#include <vector>

typedef std::vector<double> Slice;

/**
 * Three-dimensional tensor of size N1 * N2 * N3
 */
class Grid3D {
public:
    int size;
    int shape[3];

    explicit Grid3D(int N);

    Grid3D(int N1, int N2, int N3);

    ~Grid3D();

    double &operator()(int i, int j, int k);

    double operator()(int i, int j, int k) const;

    int getSliceSize(int axis) const;

    Slice getSlice(int index, int axis);

    void setSlice(int index, int axis, const Slice &slice);

    void writeToFile(std::ofstream &outFile) const;

private:
    int _cfI, _cfJ; // Coefficients for indexing

    std::vector<double> raw;

    inline int _index(int i, int j, int k) const { return i * _cfI + j * _cfJ + k; }

    void init();
};