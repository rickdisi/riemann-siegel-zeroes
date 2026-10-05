#pragma once

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "correction_terms.hpp"

// Riemann-Siegel theta function: arg Gamma(1/4 + it/2) - (t/2) ln(pi)
// Computed via its Stirling expansion, valid for t >> 1:
// theta(t) = (t/2) ln(t/2pi) - t/2 - pi/8 + 1/(48t) + 7/(5760 t^3) + ...
inline double theta(double t) {

    double result = (t / 2.0)
                    * std::log(t / (2.0 * pi))
                    - (t / 2.0)
                    - (pi / 8.0)
                    + (1.0 / (48.0 * t))
                    + (7.0 / (5760.0 * t * t * t));

    return result;
}

// Main sum of the Riemann-Siegel formula
// N = floor(sqrt(t/2pi)) is the balance point of the approximate functional equation.
// zMain(t) = 2 * sum_{n=1}^{N} cos(theta(t) - t*ln(n)) / sqrt(n)
inline double zMain(double t) {

    long N = std::floor(std::sqrt(t / (2.0 * pi))); // std::floor is technically redundant but kept for clarity

    double sum = 0.0;
    for (int n = 1; n <= N; ++n) {
        double term = std::cos(theta(t) - t * std::log(n)) / std::sqrt(n); // Summation term
        sum += term;
    }
    double mainTerm = 2.0 * sum;

    return mainTerm;
}

// Full Z(t): main sum plus the correction terms.
// Z(t) = zMain(t) + (-1)^(N-1) * (2*pi/t)^(1/4) * [c0(p) + c1(p)*(2*pi/t)^(1/2) + c2(p)*(2*pi/t)]
// N and p are computed the same way as in zMain.
inline double z(double t) {

    long N = std::floor(std::sqrt(t / (2.0 * pi))); // std::floor is technically redundant but kept for clarity
    double p = std::sqrt(t / (2.0 * pi)) - N;

    double sign = ((N - 1) % 2 == 0) ? 1.0 : -1.0; // Shorthand for if (N-1) even, then sign == +1, if odd then -1.

    double result =
        zMain(t)
        + sign * (std::pow((2.0 * pi) / t, 1.0 / 4.0)) // replaced std::pow(-1, N-1) with cheaper "sign" variable.
        * (
            c0(p) 
            + (c1(p) * std::sqrt(2.0 * pi / t))
            + c2(p) * (2.0 * pi / t)
        );

    return result;
}

// Scans z(t) from tMin to tMax in steps of `step` and returns every bracket
// (t, t+step) where z(t) changes sign.
inline std::vector<std::pair<double, double>> findSignChanges(double tMin, double tMax, double step) {
    
    std::vector<std::pair<double, double>> brackets;

    double tPrev = tMin;
    double zPrev = z(tPrev);

    while (tPrev < tMax) {
        double tNew = tPrev + step;
        double zNew = z(tNew);

        if (zPrev * zNew < 0.0) {
            brackets.push_back({tPrev, tNew});
        }

        tPrev = tNew;
        zPrev = zNew;
    }
    
    return brackets;
}

// Refines a bracket [a, b] (where z(a) and z(b) have opposite signs) down to
// a precise root via bisection. Each iteration halves the bracket.
inline double bisect(double a, double b, int iterations) {

    double za = z(a);

    for (int i = 0; i < iterations; ++i) {
        double m = (a + b) / 2.0;
        double zm = z(m);

        if (za * zm < 0.0) {
            b = m;
        }
        else {
            a = m;
            za = zm;
        }
        
    }
    return (a + b) / 2.0;
}

// Estimates a tMax large enough to contain at least `targetCount` zeros,
// using the Riemann-von Mangoldt formula N(T) ~ (T/2pi)*ln(T/2pi*e) + 7/8
// as an average-case estimate.
inline double estimateTMax(int targetCount) {

    auto N = [](double T) {
        double nOfT = (T / (2.0 * pi)) * std::log(T/(2.0 * pi * std::exp(1.0))) + (7.0 / 8.0);
        
        return nOfT;
    };

    double T = 100.0;
    while (N(T) < targetCount) {
        T = T * 2.0;
    }

    T = T * 1.1;

    return T;
}

// Finds the n-th Gram point g_n: the t where theta(t) = n * pi
// Scoped to n >= 0.
inline double gram(int n) {
    double target = n * pi;

    // Stage 1: doubling search for a bracket, same shape as estimateTMax
    double t = 10.0;
    while (theta(t) < target) {
        t = t * 2.0;
    }
    double a = t / 2.0; // the previous t, before the last doubling
    double b = t; // the t that first satisfied the condition

    // Stage 2: bisection, same shape as bisect(). Comparing theta(x) - target instead of z(x)
    double fa = theta(a) - target;

    for (int i = 0; i < 50; ++i) {
        double m = (a + b) / 2.0;
        double fm = theta(m) - target;

        if (fa * fm < 0.0) {
            b = m;
        } else {
            a = m;
            fa = fm;
        }
    }

    return (a + b) / 2.0;
}

// Checks Gram's Law at n: whether (-1)^n * Z(g_n) > 0, i.e. z(gram(n))
// alternates sign as n increases by one. True at most n.
// A false return needs Gram-block handling rather than the one zero per interval assumption
inline bool satisfiesGramLaw(int n) {

    double zg = z(gram(n));
    double sign = (n % 2 == 0) ? 1.0 : -1.0; // Shorthand for if n even, then sign == +1, if odd then -1.

    if (sign * zg > 0.0) {
        return true;
    }

    return false;
}

// Verifies a Gram block [gram(nStart), gram(nEnd)]: checks that exactly
// (nEnd - nStart) sign changes of z(t) occur across the whole span, which
// is the completeness condition for a block (individual sub-intervals
// don't need exactly one each, just the block total).
//
// TODO(Riccardo): implement this.
inline bool verifyGramBlock(int nStart, int nEnd, double step) {
    return false;
}
