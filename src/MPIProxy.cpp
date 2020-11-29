//
// Created by Dmitry Lyukov on 29.11.2020.
//

#include "MPIProxy.h"
#include "utils.h"
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

void MPIProxy::sendDouble(double value, int receiver) const {
    MPI_Send(&value, 1, MPI_DOUBLE,receiver, 0, MPI_COMM_WORLD);
}

double MPIProxy::receiveDouble(int sender) const {
    MPI_Status dummyStatus;
    double value;
    MPI_Recv(&value, 1, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD, &dummyStatus);
    return value;
}

double MPIProxy::maxOverAll(double value) const {
    double result;
    MPI_Reduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, getMainProcId(), MPI_COMM_WORLD);
    return result;
}

int MPIProxy::getNumOfProcessors() const {
    int n_proc;
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    return n_proc;
}

void MPIProxy::barrier() const {
    MPI_Barrier(MPI_COMM_WORLD);
}

bool MPIProxy::isMainProcess() const {
    return getRank() == getMainProcId();
}

double MPIProxy::time() const {
    return MPI_Wtime();
}
