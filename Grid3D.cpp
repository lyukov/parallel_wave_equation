#include "Grid3D.h"

#include "log.h"

Grid3D::Grid3D(int N1, int N2, int N3) : size(N1 * N2 * N3), _cfI(N2 * N3), _cfJ(N3) {
    shape[0] = N1;
    shape[1] = N2;
    shape[2] = N3;
    init();
}

Grid3D::Grid3D(int N) : size(N * N * N), _cfI(N * N), _cfJ(N) {
    shape[0] = N;
    shape[1] = N;
    shape[2] = N;
    init();
}

void Grid3D::init() {
    raw = new double[size];
    LOG << "Created raw array of size " << size << std::endl;
}

Grid3D::~Grid3D() {
    if (raw) {
        delete[] raw;
        LOG << "Deleted raw array of size " << size << std::endl;
    }
}

double &Grid3D::operator()(int i, int j, int k) {
    return raw[_index(i, j, k)];
}

double Grid3D::operator()(int i, int j, int k) const {
    return raw[_index(i, j, k)];
}

void Grid3D::writeToFile(std::ofstream &outFile) const {
    outFile.write((char *) raw, size * sizeof(double));
}

int Grid3D::getSliceSize(int axis) const {
    return size / shape[axis];
}

Slice Grid3D::getSlice(int index, int axis) {
    Slice slice = Slice(getSliceSize(axis));
    int idx = 0;
    if (axis == 0) {
        for (int j = 0; j < shape[1]; ++j) {
            for (int k = 0; k < shape[2]; ++k) {
                slice[idx++] = (*this)(index, j, k);
            }
        }
    } else if (axis == 1) {
        for (int i = 0; i < shape[0]; ++i) {
            for (int k = 0; k < shape[3]; ++k) {
                slice[idx++] = (*this)(i, index, k);
            }
        }
    } else {
        for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                slice[idx++] = (*this)(i, j, index);
            }
        }
    }
    return slice;
}

void Grid3D::setSlice(int index, int axis, const Slice &slice) {
    int idx = 0;
    if (axis == 0) {
        for (int j = 0; j < shape[1]; ++j) {
            for (int k = 0; k < shape[2]; ++k) {
                (*this)(index, j, k) = slice[idx++];
            }
        }
    } else if (axis == 1) {
        for (int i = 0; i < shape[0]; ++i) {
            for (int k = 0; k < shape[2]; ++k) {
                (*this)(i, index, k) = slice[idx++];
            }
        }
    } else {
        for (int i = 0; i < shape[0]; ++i) {
            for (int j = 0; j < shape[1]; ++j) {
                (*this)(i, j, index) = slice[idx++];
            }
        }
    }
}
