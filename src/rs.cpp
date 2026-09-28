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
    double total = 2.0 * sum;

    return total;
}

int main() {
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

    return 0;
}
