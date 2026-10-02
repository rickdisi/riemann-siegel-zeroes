#include "rs.hpp"

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

    // Verification: refine the first 10 brackets with bisect() and compare
    // against the published first 10 zeros, to 8 decimal places.
    double knownZeros[10] = {
        14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588,
        37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478
    };
    std::printf("\nrefined zeros vs known values:\n");
    for (int i = 0; i < 10; ++i) {
        double refined = bisect(brackets[i].first, brackets[i].second, 50);
        double error = refined - knownZeros[i];
        std::printf("  zero %d: refined=%.9f  known=%.9f  error=%.2e\n", i + 1, refined, knownZeros[i], error);
    }

    return 0;
}
