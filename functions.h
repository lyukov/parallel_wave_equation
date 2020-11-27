#pragma once

class Function4D {
   public:
    virtual double operator()(double t, double x, double y, double z) const = 0;
};

class Function3D {
   public:
    virtual double operator()(double x, double y, double z) const = 0;
};