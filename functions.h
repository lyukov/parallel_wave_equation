#pragma once

class Function4D {
public:
    virtual double operator()(double t, double x, double y, double z) const = 0;
};

class Function3D {
public:
    virtual double operator()(double x, double y, double z) const = 0;
};

class U : public Function4D {
    const double L_x;
    const double L_y;
    const double L_z;
    const double a_t;

public:
    U(double L_x, double L_y, double L_z);

    double operator()(double t, double x, double y, double z) const;
};

class Phi : public Function3D {
    double L_x;
    double L_y;
    double L_z;

public:
    Phi(double L_x, double L_y, double L_z);

    double operator()(double x, double y, double z) const;
};