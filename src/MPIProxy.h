#pragma once

#include <vector>

class MPIProxy {
public:
    MPIProxy(int *argc, char ***argv);

    ~MPIProxy();

    int getRank() const;

    int getNumOfProcessors() const;

    void sendVector(const std::vector<double> &data, int receiver) const;

    std::vector<double> receiveVector(int size, int sender) const;

    double maxOverAll(double value) const;

    int getMainProcId() const { return 0; }

private:
    MPIProxy(const MPIProxy &); // non construction-copyable
    MPIProxy &operator=(const MPIProxy &); // non copyable
};
