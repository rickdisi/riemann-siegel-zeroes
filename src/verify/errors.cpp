#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>

#include "rs.hpp"

// Z(t) keeping the first `terms` correction terms: C0 only (1), C0..C1 (2), C0..C2 (3).
// Same formula and evaluation order as z(t) in rs.hpp, which is left untouched;
// main() below checks that zTerms(t, 3) agrees with z(t) before anything is written.
inline double zTerms(double t, int terms) {

    long N = std::floor(std::sqrt(t / (2.0 * pi)));
    double p = std::sqrt(t / (2.0 * pi)) - N;

    double sign = ((N - 1) % 2 == 0) ? 1.0 : -1.0;

    double bracket = c0(p);
    if (terms >= 2) {
        bracket += c1(p) * std::sqrt(2.0 * pi / t);
    }
    if (terms >= 3) {
        bracket += c2(p) * (2.0 * pi / t);
    }

    return zMain(t) + sign * std::pow((2.0 * pi) / t, 1.0 / 4.0) * bracket;
}

int main() {

    // Heights are set by N0 = floor(sqrt(t / 2pi)), so t ~ 2*pi*N0^2 runs from 25 to ~2.5e5.
    // Above ~3e5, double-precision rounding of the phase (~1e-9) approaches the truncation error.
    const int n0List[] = {2, 3, 5, 8, 12, 20, 35, 60, 100, 150, 200};
    const int pCount = 100;

    // Sample p across [0, 1) directly: stepping in t barely moves p at large t.
    // Points at p = (i + 0.5) / 100 keep clear of the c1/c2 singularities at p = 0.25, 0.75.
    double maxDiff = 0.0;
    for (int n0 : n0List) {
        for (int i = 0; i < pCount; ++i) {
            double t = 2.0 * pi * std::pow(n0 + (i + 0.5) / pCount, 2.0);
            maxDiff = std::max(maxDiff, std::fabs(zTerms(t, 3) - z(t)));
        }
    }
    std::printf("max |zTerms(t, 3) - z(t)| over all sample points: %.3e\n", maxDiff);
    if (maxDiff > 1e-12) {
        std::fprintf(stderr, "FAILED: zTerms(t, 3) does not match z(t), so the experiment would not measure the production code\n");
        return 1;
    }

    std::ofstream out("data/term_errors.csv");
    if (!out) {
        std::fprintf(stderr, "failed to open data/term_errors.csv for writing\n");
        return 1;
    }

    out << "N,p,t,z1,z2,z3\n";
    out << std::setprecision(17); // the default 6 digits would hide the errors being measured

    for (int n0 : n0List) {
        for (int i = 0; i < pCount; ++i) {
            double p = (i + 0.5) / pCount;
            double t = 2.0 * pi * std::pow(n0 + p, 2.0);
            out << n0 << "," << p << "," << t << ","
                << zTerms(t, 1) << "," << zTerms(t, 2) << "," << zTerms(t, 3) << "\n";
        }
    }

    std::printf("wrote data/term_errors.csv (%zu heights x %d values of p)\n", sizeof(n0List) / sizeof(n0List[0]), pCount);
    return 0;
}
