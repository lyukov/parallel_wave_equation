#include <stdio.h>

typedef double(Function)(double t, double x, double y, double z);

class Grid {
    double T;
    double L_x;
    double L_y;
    double L_z;
    Function u;

   public:
    Grid(double T, double L_x, double L_y, double L_z, Function u) {}
};

int main(int argc, char **argv) {
    return 0;
}