#include "Grid3D.h"

#include "utils.h"

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
    raw.resize(size);
}

double &Grid3D::operator()(int i, int j, int k) {
    return raw[_index(i, j, k)];
}

double Grid3D::operator()(int i, int j, int k) const {
    return raw[_index(i, j, k)];
}

void Grid3D::writeToFile(const std::string &fileName) const {
    std::ofstream outFile(fileName.c_str(), std::ios::out | std::ios::binary);
    outFile.write((char *) raw.data(), size * sizeof(double));
    outFile.close();
}

int Grid3D::getSliceSize(int axis) const {
    return size / shape[axis];
}

void Grid3D::getSliceParams(int axis, int &c0, int &c1, int &c2, int &N, int &M) const {
    int cf[3] = {shape[1] * shape[2], shape[2], 1};
    c0 = cf[axis % 3];
    c1 = cf[(axis + 1) % 3];
    c2 = cf[(axis + 2) % 3];
    N = shape[(axis + 1) % 3];
    M = shape[(axis + 2) % 3];
}

std::vector<double> Grid3D::getSlice(int index, int axis) {
    std::vector<double> slice = std::vector<double>(getSliceSize(axis));
    int c0, c1, c2, N, M;
    getSliceParams(axis, c0, c1, c2, N, M);
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        int idx = c0 * index + c1 * i;
        int slice_idx = i * M;
        for (int j = 0; j < M; ++j) {
            slice[slice_idx] = raw[idx];
            ++slice_idx;
            idx += c2;
        }
    }
    return slice;
}

void Grid3D::setSlice(int index, int axis, const std::vector<double> &slice) {
    int c0, c1, c2, N, M;
    getSliceParams(axis, c0, c1, c2, N, M);
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        int idx = c0 * index + c1 * i;
        int slice_idx = i * M;
        for (int j = 0; j < M; ++j) {
            raw[idx] = slice[slice_idx];
            ++slice_idx;
            idx += c2;
        }
    }
}

void Grid3D::setZeros(int index, int axis) {
    int c0, c1, c2, N, M;
    getSliceParams(axis, c0, c1, c2, N, M);
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        int idx = c0 * index + c1 * i;
        for (int j = 0; j < M; ++j) {
            raw[idx] = 0;
            idx += c2;
        }
    }
}

std::vector<double> &Grid3D::getFlatten() {
    return raw;
}
