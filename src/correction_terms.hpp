#pragma once

#include <cmath>
#include <vector>

constexpr double pi = 3.14159265358979323846;

// Taylor coefficients of C1(z) expanded around z = -0.5.
const std::vector<double> c1TaylorCoeffs = {
    - 0.010416666667,
    0.003510473973,
    0.074180889406,
    - 0.126319511253,
    0.122624182326,
    - 0.053783813433,
    - 0.020927470267,
    0.048901792109,
    - 0.042419263336,
    0.018941888162,
};

// Taylor coefficients of C2(z) expanded around z = -0.5.
const std::vector<double> c2TaylorCoeffs = {
    0.004612789401,
    0.004611927296,
    - 0.012307171355,
    0.007943242061,
    0.019669481904,
    - 0.041777952238,
    0.033580528430,
    - 0.009131134250,
    - 0.012724251032,
    0.018527514859,
};

// Evaluates a Taylor series sum_i coeffs[i] * (x - center)^i via Horner's
// method (processes highest-degree term first, avoiding repeated pow() calls).
inline double evalTaylorSeries(const std::vector<double>& coeffs, double center, double x) {

    double dx = x - center;
    double result = coeffs.back();

    for (int i = coeffs.size() - 2; i >= 0; --i) {
        result = result * dx + coeffs[i];
    }

    return result;
}

// First Riemann-Siegel correction term.
// p is the fractional leftover of N: p = sqrt(t/2pi) - N
// c0(p) = cos(2*pi*(p^2 - p - 1/16)) / cos(2*pi*p)
inline double c0(double p) {

    double c0Coeff = std::cos(2.0 * pi * (p * p - p - (1.0 / 16.0))) / std::cos(2.0 * pi * p);

    return c0Coeff;
}

// Closed form evaluation of C1(z): F'''(z) / (12*pi^2).
// p is the same fractional leftover as in c0. Unstable near z = 0.5
// Use the Taylor expansion (c1TaylorCoeffs) close to z = ±0.5.
inline double c1Direct(double p) {

    double zSubst = 1.0 - 2.0 * p;
    double g = (pi / 2.0) * (zSubst * zSubst + 3.0 / 4.0);
    double sinPiZ = std::sin(pi * zSubst);
    double cosPiZ = std::cos(pi * zSubst);
    double sinG = std::sin(g);
    double cosG = std::cos(g);

    double zSubst2 = zSubst * zSubst;
    double sinPiZ2 = sinPiZ * sinPiZ;
    double cosPiZ2 = cosPiZ * cosPiZ;
    double cosPiZ3 = cosPiZ2 * cosPiZ;
    double cosPiZ4 = cosPiZ2 * cosPiZ2;

    double term1 = zSubst * (pi * zSubst2 * sinG - 3.0 * cosG) * cosPiZ3;
    double term2 = -3.0 * pi * zSubst * (sinPiZ2 + 1.0) * sinG * cosPiZ;
    double term3 = -3.0 * (pi * zSubst2 * cosG + sinG) * sinPiZ * cosPiZ2;
    double term4 = pi * (sinPiZ2 + 5.0) * sinPiZ * cosG;

    double c1Coeff = (term1 + term2 + term3 + term4) / (12.0 * cosPiZ4);

    return c1Coeff;
}

// First/second Riemann-Siegel correction term (C1). p is the same
// fractional leftover as in c0. Dispatches to the Taylor expansion
// (c1TaylorCoeffs) within 0.001 of the singular points z = ±0.5, and to
// the direct closed form (c1Direct) otherwise.
inline double c1(double p) {

    double zSubst = 1.0 - 2.0 * p;
    double result;

    if (- 0.5 - 0.001 <= zSubst && zSubst <= - 0.5 + 0.001) {
        result = evalTaylorSeries(c1TaylorCoeffs, -0.5, zSubst);
    }
    else if (0.5 - 0.001 <= zSubst && zSubst <= 0.5 + 0.001) {
        result = - evalTaylorSeries(c1TaylorCoeffs, - 0.5, -zSubst);
    }
    else {
        result = c1Direct(p);
    }

    return result;
}

// Closed form evaluation of C2(z): d_2,0 * F^(6)(z)/pi^4 + d_2,1 * F''(z)/pi^2.
// Unstable near z = ±0.5. Use c2Taylor close to z = ±0.5.
inline double c2Direct(double p) {

    double zSubst = 1.0 - 2.0 * p;
    double g = (pi / 2.0) * (zSubst * zSubst + 3.0 / 4.0);
    double sinPiZ = std::sin(pi * zSubst);
    double cosPiZ = std::cos(pi * zSubst);
    double sinG = std::sin(g);
    double cosG = std::cos(g);

    double pi2 = pi * pi;
    double pi3 = pi2 * pi;
    double cosPiZ2 = cosPiZ * cosPiZ;
    double cosPiZ3 = cosPiZ2 * cosPiZ;
    double cosPiZ4 = cosPiZ2 * cosPiZ2;
    double cosPiZ5 = cosPiZ4 * cosPiZ;
    double cosPiZ6 = cosPiZ4 * cosPiZ2;
    double cosPiZ7 = cosPiZ6 * cosPiZ;
    double sinPiZ2 = sinPiZ * sinPiZ;
    double sinPiZ4 = sinPiZ2 * sinPiZ2;
    double sinPiZ6 = sinPiZ4 * sinPiZ2;
    double zSubst2 = zSubst * zSubst;
    double zSubst4 = zSubst2 * zSubst2;
    double zSubst6 = zSubst4 * zSubst2;

    double term0 = pi3;
    double term1 = sinPiZ4;
    double term2 = sinPiZ2;
    double term3 = sinPiZ * zSubst;
    double term4 = cosPiZ * term3;
    double term5 = pi2;
    double term6 = 3.0 * cosG;
    double term7 = zSubst2;
    double term8 = cosPiZ6;
    double term9 = cosG * term0;
    double term10 = pi * cosG * term7;
    double term11 = cosPiZ2 * (sinG + term10);
    double term12 = 15.0 * sinG;
    double term13 = term5 * zSubst4;
    double term14 = cosPiZ4;
    double term15 = term2 + 1.0;
    double term16 = pi * sinG;

    double c2Coeff = (1.0 / 288.0) * (
        6.0 * pi * cosPiZ5 * term3 * (-sinG * term13 + 10.0 * term10 + term12)
        + 20.0 * cosPiZ3 * term3 * term5 * (term2 + 5.0) * (pi * sinG * term7 - term6)
        - 6.0 * sinG * term0 * term4 * (term1 + 58.0 * term2 + 61.0)
        - 15.0 * term11 * term5 * (term1 + 18.0 * term2 + 5.0)
        + 15.0 * pi * term14 * term15 * (cosG * term13 + 6.0 * term16 * term7 - term6)
        + 18.0 * term14 * (pi * cosG * term15 - term11 - 2.0 * term16 * term4)
        + term8 * (45.0 * pi * cosG * term7 + 15.0 * sinG - term12 * term13 - term9 * zSubst6)
        + term9 * (62.0 * sinPiZ6 - 4.0 * term1 + 662.0 * term2 + 61.0 * term8)
    ) / (pi * cosPiZ7);

    return c2Coeff;
}

// Third Riemann-Siegel correction term. p is the same fractional leftover
// as in c0/c1. Dispatches to the Taylor expansion (c2TaylorCoeffs) within
// 0.02 of the singular points z=±0.5, and to the direct closed form c2Direct otherwise.
inline double c2(double p) {

    double zSubst = 1.0 - 2.0 * p;
    double result;

    if (- 0.5 - 0.02 <= zSubst && zSubst <= - 0.5 + 0.02) {
        result = evalTaylorSeries(c2TaylorCoeffs, -0.5, zSubst);
    }
    else if (0.5 - 0.02 <= zSubst && zSubst <= 0.5 + 0.02) {
        result = evalTaylorSeries(c2TaylorCoeffs, -0.5, -zSubst);
    }
    else {
        result = c2Direct(p);
    }

    return result;
}