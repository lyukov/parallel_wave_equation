//
// Created by Dmitry Lyukov on 29.11.2020.
//

#include "MPIProxy.h"

#include <mpi.h>

MPIProxy::MPIProxy(int *argc, char ***argv) {
    MPI_Init(argc, argv);
}

MPIProxy::~MPIProxy() {
    MPI_Finalize();
}

int MPIProxy::getRank() const {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

void MPIProxy::sendVector(const std::vector<double> &data, int receiver) const {
    MPI_Request dummyRequest;
    MPI_Isend(data.data(), data.size(), MPI_DOUBLE, receiver, 0, MPI_COMM_WORLD, &dummyRequest);
}

std::vector<double> MPIProxy::receiveVector(int size, int sender) const {
    MPI_Status dummyStatus;
    std::vector<double> data(size);
    MPI_Recv(data.data(), size, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD, &dummyStatus);
    return data;
}

double MPIProxy::maxOverAll(double value) const {
    double result;
    MPI_Reduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return result;
}
