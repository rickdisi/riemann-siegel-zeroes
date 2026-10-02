#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>

constexpr double pi = 3.14159265358979323846;

// Riemann-Siegel theta function: arg Gamma(1/4 + it/2) - (t/2) ln(pi)
// Computed via its Stirling expansion, valid for t >> 1:
// theta(t) = (t/2) ln(t/2pi) - t/2 - pi/8 + 1/(48t) + 7/(5760 t^3) + ...
double theta(double t) {

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
double zMain(double t) {

    long N = std::floor(std::sqrt(t / (2.0 * pi))); // std::floor is technically redundant but kept for clarity

    double sum = 0.0;
    for (int n = 1; n <= N; ++n) {
        double term = std::cos(theta(t) - t * std::log(n)) / std::sqrt(n); // Summation term
        sum += term;
    }
    double mainTerm = 2.0 * sum;

    return mainTerm;
}

// First Riemann-Siegel correction term.
// p is the fractional leftover of N: p = sqrt(t/2pi) - N
// c0(p) = cos(2*pi*(p^2 - p - 1/16)) / cos(2*pi*p)
double c0(double p) {

    double c0Coeff = std::cos(2.0 * pi * (p * p - p - (1.0 / 16.0))) / std::cos(2.0 * pi * p);

    return c0Coeff;
}

// Second Riemann-Siegel correction term.
// p is the same fractional leftover as in c0.
// C1 is the third derivative of F(z), scaled by a constant: C1(z) = F'''(z) / (12*pi^2).
double c1(double p) {

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

// Third Riemann-Siegel correction term.
// p is the same fractional leftover as in c0/c1.
// C2(z) = d_2,0 * F^(6)(z)/pi^4 + d_2,1 * F''(z)/pi^2.
double c2(double p) {

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

// Full Z(t): main sum plus the correction terms.
// Z(t) = zMain(t) + (-1)^(N-1) * (2*pi/t)^(1/4) * [c0(p) + c1(p)*(2*pi/t)^(1/2) + c2(p)*(2*pi/t)]
// N and p are computed the same way as in zMain.
double z(double t) {

    long N = std::floor(std::sqrt(t / (2.0 * pi))); // std::floor is technically redundant but kept for clarity
    double p = std::sqrt(t / (2.0 * pi)) - N;

    double sign = ((N - 1) % 2 == 0) ? 1.0 : -1.0; // Shorthand for if (N-1) even, then sign==+1, if odd then -1.

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
std::vector<std::pair<double, double>> findSignChanges(double tMin, double tMax, double step) {
    
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

int main() {
    // Tests written by Claude
    // Sanity check: g0 = 17.8455995405 is the first Gram point, defined by
    // theta(g0) = 0 (confirmed against Wikipedia's Gram point article).
    double g0 = 17.8455995405;
    std::printf("theta(g0) = %.10f  (expected close to 0)\n", theta(g0));

    // Sanity check: zMain should bracket the first known zero of zeta on the
    // critical line, t = 14.134725... Correction terms are not added yet, and
    // at this height N = floor(sqrt(t/2pi)) = 1, so the main sum alone is a
    // crude approximation: it gets a sign change in roughly the right place,
    // but not the precise location (confirmed by scanning zMain t=10..55:
    // it crosses zero near all 10 known zeroes, shifted by up to ~1 in t).
    // Widened bracket to actually catch the sign change at this precision.
    double before = 14.0;
    double after = 15.0;
    double zBefore = zMain(before);
    double zAfter = zMain(after);
    std::printf("zMain(%.1f) = %.6f\n", before, zBefore);
    std::printf("zMain(%.1f) = %.6f\n", after, zAfter);
    if (zBefore * zAfter < 0.0) {
        std::printf("sign change detected -> a zero lies between %.1f and %.1f\n", before, after);
    } else {
        std::printf("no sign change -> check zMain\n");
    }

    // Sanity check: z(t) (main sum + c0 correction) should land much closer
    // to the true first zero, t1 = 14.134725142..., than zMain alone did.
    // zMain only guaranteed a sign change somewhere in [14.0, 15.0]; z should
    // be small right at t1 itself, and should also catch a sign change in a
    // much narrower bracket around t1.
    double t1 = 14.134725142;
    std::printf("z(t1) = %.6f  (expected much closer to 0 than zMain alone)\n", z(t1));

    double zBeforeNarrow = z(14.10);
    double zAfterNarrow = z(14.17);
    std::printf("z(14.10) = %.6f\n", zBeforeNarrow);
    std::printf("z(14.17) = %.6f\n", zAfterNarrow);
    if (zBeforeNarrow * zAfterNarrow < 0.0) {
        std::printf("sign change detected in narrow bracket -> c0 correction improved the location\n");
    } else {
        std::printf("no sign change in narrow bracket -> check c0/z\n");
    }

    // Sanity check: scanning t=10..55 with step=0.1 should find exactly 10
    // sign-change brackets, one per known zero: 14.134725, 21.022040,
    // 25.010858, 30.424876, 32.935062, 37.586178, 40.918719, 43.327073,
    // 48.005151, 49.773832.
    auto brackets = findSignChanges(10.0, 55.0, 0.1);
    std::printf("found %zu sign-change brackets:\n", brackets.size());
    for (auto& bracket : brackets) {
        std::printf("  (%.2f, %.2f)\n", bracket.first, bracket.second);
    }

    return 0;
}
