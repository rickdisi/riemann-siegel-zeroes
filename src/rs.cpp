#include <cmath>
#include <cstdio>

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
// p is the fractional leftoverof N: p = sqrt(t/2pi) - N
// c0(p) = cos(2*pi*(p^2 - p - 1/16)) / cos(2*pi*p)
double c0(double p) {

    double c0Coeff = std::cos(2.0 * pi * (p * p - p - (1.0 / 16.0))) / std::cos(2 * pi * p);

    return c0Coeff;
}

// Second Riemann-Siegel correction term.
// p is the same fractional leftover as in c0.
// C1 is the third derivative of F(z), scaled by a constant: C1(z) = F'''(z) / (12*pi^2).
double c1(double p) {

    double zSubst = 1 - 2 * p;
    double g = (pi / 2.0) * (zSubst * zSubst + 3.0/4.0);

    double c1Coeff = ( 
        zSubst * (pi * zSubst * zSubst * std::sin(g) - 3 * std::cos(g)) 
        * std::cos(pi * zSubst) * std::cos(pi * zSubst) * std::cos(pi * zSubst)
        - 3 * pi * zSubst * (std::sin(pi * zSubst) * std::sin(pi * zSubst) + 1)
        * std::sin(g) * std::cos(pi * zSubst)
        - 3 * (pi * zSubst * zSubst * std::cos(g) + std::sin(g)) 
        * std::sin(pi * zSubst) * std::cos(pi * zSubst) * std::cos(pi * zSubst)
        + pi * (std::sin(pi * zSubst) * std::sin(pi * zSubst) + 5) 
        * std::sin(pi * zSubst) * std::cos(g)
    )
        / (
            12.0 * std::cos(pi * zSubst) * std::cos(pi * zSubst) * std::cos(pi * zSubst) * std::cos(pi * zSubst)
        );
    
    return c1Coeff;
}

// Full Z(t): main sum plus the correction terms.
// Z(t) = zMain(t) + (-1)^(N-1) * (2*pi/t)^(1/4) * [c0(p) + c1(p)*(2*pi/t)^(1/2)]
// N and p are computed the same way as in zMain. Each further correction
// term (c2, c3, c4) would carry one more power of (2*pi/t)^(1/2) inside
// the brackets.
double z(double t) {

    long N = std::floor(std::sqrt(t / (2.0 * pi))); // std::floor is technically redundant but kept for clarity
    double p = std::sqrt(t / (2.0 * pi)) - N;

    double sign = ((N - 1) % 2 == 0) ? 1.0 : -1.0; // Shorthand for if (N-1) even, then sign==+1, if odd then -1.
    
    double result = 
        zMain(t) 
        + sign * (std::pow((2.0 * pi) / t, 1.0/4.0)) // replaced std::pow(-1, N-1) with cheaper "sign" variable.
        * (c0(p) + (c1(p) * std::sqrt(2.0 * pi / t))); 

    return result;
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

    return 0;
}
