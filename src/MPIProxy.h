#pragma once

#include <vector>

class MPIProxy {
public:
    MPIProxy(int *argc, char ***argv);

    ~MPIProxy();

    int getRank() const { return _rank; }

    int getNumOfProcessors() const { return _n_proc; }

    std::vector<int> createDims(int nnodes, int ndims) const;

    void sendVector(std::vector<double> &data, int receiver) const;

    std::vector<double> receiveVector(int size, int sender) const;

    void sendDouble(double value, int receiver) const;

    double receiveDouble(int sender) const;

    double maxOverAll(double value) const;

    double sumOverAll(double value) const;

    int getMainProcId() const { return 0; }

    bool isMainProcess() const;

    double time() const;

    void barrier() const;

private:
    int _rank;
    int _n_proc;

    MPIProxy(const MPIProxy &); // non construction-copyable
    MPIProxy &operator=(const MPIProxy &); // non copyable
};
