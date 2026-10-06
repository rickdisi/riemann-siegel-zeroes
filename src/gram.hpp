#pragma once

#include <vector>
#include <utility>

#include "rs.hpp"

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
// is the completeness condition for a block.
inline bool verifyGramBlock(int nStart, int nEnd, double step) {

    double tLower = gram(nStart);
    double tUpper = gram(nEnd);

    auto brackets = findSignChanges(tLower, tUpper, step);
    
    if (brackets.size() == static_cast<size_t>(nEnd - nStart)) { // Avoids (nEnd - nStart) negative, even if impossible.
        return true;
    }

    return false;
}

// Scans n = nMin + 1 to nMax, grouping Gram points into blocks. Each block
// is bounded by two Gram points that satisfy Gram's Law, with anything in
// between (if present) violating it.
inline std::vector<std::pair<int, int>> findGramBlocks(int nMin, int nMax) {

    int blockStart = nMin;
    std::vector<std::pair<int, int>> results;

    for (int n = nMin + 1; n <= nMax; ++n) {
        if (satisfiesGramLaw(n) == true) {
            results.push_back({blockStart, n});
            blockStart = n;
        }
    }

    return results;
}

// Verifies every Gram block in [nMin, nMax], returning the list of blocks that failed verification.
// Empty result means every block was confirmed complete.
inline std::vector<std::pair<int, int>> verifyAllGramBlocks(int nMin, int nMax, double step) {

    auto blocks = findGramBlocks(nMin, nMax);
    std::vector<std::pair<int, int>> failedBlocks;

    for (auto& block : blocks) { // Range-based loop.
        bool ok = verifyGramBlock(block.first, block.second, step);

        if (!ok) {
            failedBlocks.push_back(block);
        }
    }

    return failedBlocks;
}