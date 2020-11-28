#pragma once

class Function4D {
public:
    virtual double operator()(double t, double x, double y, double z) const = 0;
};

class Function3D {
public:
    virtual double operator()(double x, double y, double z) const = 0;
};

class U {
    const double L_x;
    const double L_y;
    const double L_z;
    const double a_t;
public:
    inline U(double L_x, double L_y, double L_z) : L_x(L_x), L_y(L_y), L_z(L_z),
            a_t(M_PI * sqrt(4.0 / (L_x * L_x) + 1.0 / (L_y * L_y) + 1.0 / (L_z * L_z))) {}

    inline double operator()(double t, double x, double y, double z) const {
        return sin(M_2_PI * x / L_x) * sin(M_PI * y / L_y) * sin(M_2_PI * z / L_z) * cos(a_t * t);
    }
};

class Phi {
    double L_x;
    double L_y;
    double L_z;
public:
    inline Phi(double L_x, double L_y, double L_z) : L_x(L_x), L_y(L_y), L_z(L_z) {}

    inline double operator()(double x, double y, double z) const {
        return sin(M_2_PI * x / L_x) * sin(M_PI * y / L_y) * sin(M_2_PI * z / L_z);
    }
};