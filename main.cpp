#include <stdio.h>

typedef double (*Function)(double t, double x, double y, double z);

class Grid {
    /* Parameters of continuous problem */
    double T;
    double L_x;
    double L_y;
    double L_z;

    /* Ground truth solution */
    Function u;

    /* Grid size is K * N * N * N */
    int N;
    int K;

    /* Grid steps */
    double h_x;
    double h_y;
    double h_z;
    double tau;

   public:
    Grid(double T, double L_x, double L_y, double L_z, Function u, int N, int K) {
        this->T = T;
        this->L_x = L_x;
        this->L_y = L_y;
        this->L_z = L_z;
        this->u = u;

        h_x = L_x / N;
        h_y = L_y / N;
        h_z = L_z / N;
        tau = T / K;
    }
};

int main(int argc, char **argv) {
    return 0;
}