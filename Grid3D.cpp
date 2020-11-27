#include "Grid3D.h"
#include "log.h"

Grid3D::Grid3D(int N) : _N1(N + 1), _N2(N + 1), _N3(N + 1), size(_N1 * _N2 * _N3) {
    raw = new double[size];
    LOG << "Created raw array of size " << size << std::endl;
}

Grid3D::Grid3D(int N1, int N2, int N3)
    : _N1(N1 + 1), _N2(N2 + 1), _N3(N3 + 1), size(_N1 * _N2 * _N3) {
    raw = new double[size];
    LOG << "Created raw array of size " << size << std::endl;
}

Grid3D::~Grid3D() {
    if (raw) {
        delete raw;
        LOG << "Deleted raw array of size " << size << std::endl;
    }
}

double &Grid3D::operator()(int i, int j, int k) {
    int idx = i * (_N2 * _N3) + j * _N3 + k;
    return raw[idx];
}

double Grid3D::operator()(int i, int j, int k) const {
    int idx = i * (_N2 * _N3) + j * _N3 + k;
    return raw[idx];
}

void Grid3D::writeToFile(std::ofstream &outFile) const {
    outFile.write((char *)raw, size * sizeof(double));
}