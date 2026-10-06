#pragma once

#include <cmath>
#include <vector>
#include <stdexcept>

#include "gram.hpp"

const double turingTMin = 168.0 * pi;

// Finds h such that (-1)^j * z(gram(j) + h) > 0, for a bad Gram point j.
// Returns 0.0 if g_j is already good.
inline double findGramOffset(int j, double step) {

    if (satisfiesGramLaw(j)) {
        return 0.0;
    }
    
    double sign = (j % 2 == 0) ? 1.0 : -1.0; // Shorthand for if j even, then sign == +1, if odd then -1.
    double gj = gram(j);
    const int maxTries = static_cast<int>(1.0 / step);

    for (int i = 1; i <= maxTries; ++i) {
        double h = i * step;

        if (sign * z(gj + h) > 0) {
            return h;
        }
        else if (sign * z(gj - h) > 0) {
            return -h;
        }
    }
    
    throw std::runtime_error("passed maxTries");
}

// Upper bound on S(g_m) from Turing's method:
// 1 + (g_{m+k} - g_m)^-1 * [2.30 + 0.128 * ln(g_{m+k} / (2*pi)) + sumH]
// sumH is the signed sum of findGramOffset(j, step) for j = m+1 .. m+k-1.
inline double turingBound(int m, int k, double sumH) {

    double gStart = gram(m);
    double gEnd = gram(m + k);

    double width = gEnd - gStart;
    double bracket = 2.30 + 0.128 * std::log(gEnd / (2.0 * pi)) + sumH;

    return 1.0 + bracket / width;
}

// Tries to prove N(g_m) = m + 1 using the span [g_m, g_{m+k}].
// Returns false if any precondition fails or the bound is not below 2.
// False means "could not prove this span".
inline bool proveSpan(int m, int k, double step) {

    if (k < 1) {
        return false;
    }
    if (gram(m) <= turingTMin) { // Not valid below this number
        return false;
    }
    if (!satisfiesGramLaw(m) || !satisfiesGramLaw(m + k)) {
        return false;
    }

    double sumH = 0.0;
    double prev = gram(m); // previous corrected point, for the strictly-increasing check

    // Interior Gram points j = m+1 .. m+k-1
    for (int j = m + 1; j <= m + k - 1; ++j) {
        double h;
        try {
            h = findGramOffset(j, step);
        } 
        catch (const std::runtime_error&) {
            return false; // No valid offset found, so this span can't be proven
        }
        double corrected = gram(j) + h;

        if (corrected <= prev) { // g_j + h_j must be strictly increasing
            return false;
        }

        sumH += h;
        prev = corrected;
    }

    // The last interior point must also sit below the right endpoint
    if (gram(m + k) <= prev) {
        return false;
    }

    return turingBound(m, k, sumH) < 2.0;
}

// Finds the smallest k in [1, kMax] for which proveSpan(m, k, step) succeeds.
// Returns that k, or -1 if no k up to kMax works.
inline int findProvableSpan(int m, int kMax, double step) {

    for (int k = 1; k <= kMax; ++k) {
        // proveSpan already returns false when the endpoint g_{m+k} violates Gram's Law, so no separate skip is needed
        if (proveSpan(m, k, step)) {
            return k;
        }
    }

    return -1;
}

// Confirms that the span [g_m, g_{m+k}] contains exactly k zeros by scanning
// for sign changes. proveSpan proves N(g_m) = m+1 and N(g_{m+k}) = m+k+1.
// A proven span must contain exactly k zeros between its endpoints.
inline bool crossCheckSpan(int m, int k, double step) {
    return verifyGramBlock(m, m + k, step);
}

// Proves that the Gram-point span [g_m, g_{m+k}] contains exactly k zeros of
// zeta in the critical strip, and that the scan finds all of them on the critical line.
//   1. N(g_m) = m + 1 via findProvableSpan(m, kMax, step)
//   2. N(g_{m+k}) = m + k + 1 via findProvableSpan(m + k, kMax, step)
//   3. crossCheckSpan(m, k, step): the scan finds exactly k sign changes
inline bool proveCompleteSpan(int m, int k, int kMax, double step) {

    if (k < 1) {
        return false;
    }

    // left endpoint
    if (findProvableSpan(m, kMax, step) == -1) {
        return false;
    }

    // right endpoint
    if (findProvableSpan(m + k, kMax, step) == -1) {
        return false;
    }

    // scan
    return crossCheckSpan(m, k, step);
}