#include "Grid3D.h"

#include "log.h"

Grid3D::Grid3D(int N1, int N2, int N3, int offset_1, int offset_2, int offset_3)
        : _N1(N1 + 1), _N2(N2 + 1), _N3(N3 + 1), size(_N1 * _N2 * _N3),
          offset_1(offset_1), offset_2(offset_2), offset_3(offset_3) {
    init();
}

Grid3D::Grid3D(int N1, int N2, int N3)
        : _N1(N1 + 1), _N2(N2 + 1), _N3(N3 + 1), size(_N1 * _N2 * _N3),
          offset_1(0), offset_2(0), offset_3(0) {
    init();
}

Grid3D::Grid3D(int N)
        : _N1(N + 1), _N2(N + 1), _N3(N + 1), size(_N1 * _N2 * _N3),
          offset_1(0), offset_2(0), offset_3(0) {
    init();
}

void Grid3D::init() {
    raw = new double[size];
    LOG << "Created raw array of size " << size << std::endl;
}

Grid3D::~Grid3D() {
    if (raw) {
        delete raw;
        LOG << "Deleted raw array of size " << size << std::endl;
    }
}

int Grid3D::_index(int i, int j, int k) const {
    return (i - offset_1) * (_N2 * _N3) + (j - offset_2) * _N3 + (k - offset_3);
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