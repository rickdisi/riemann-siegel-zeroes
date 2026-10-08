#include "rs.hpp"
#include "gram.hpp"
#include "turing.hpp"

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

    // Find the first n where Gram's Law actually fails, since n=0..10 (above)
    // all satisfied it -- we need a real violation to design the Gram-block
    // handling against, rather than guessing a literature value from memory.
    std::printf("\nscanning for the first Gram's Law violation:\n");
    int n = 0;
    while (satisfiesGramLaw(n)) {
        ++n;
    }
    std::printf("  first violation at n=%d (g_n=%.6f)\n", n, gram(n));

    // findGramOffset: for the first few bad Gram points above 168*pi (where
    // Turing's bound is valid), sign * z(g_j + h) must be positive. A good
    // Gram point must return exactly 0.0.
    std::printf("\nfindGramOffset on bad Gram points above 168*pi:\n");
    int shown = 0;
    for (int j = 290; j < 600 && shown < 4; ++j) {
        if (!satisfiesGramLaw(j)) {
            double h = findGramOffset(j, 0.01);
            double sign = (j % 2 == 0) ? 1.0 : -1.0;
            std::printf("  n=%d g_n=%.4f h=%+.2f sign*z(g_n+h)=%.4f\n", j, gram(j), h, sign * z(gram(j) + h));
            ++shown;
        }
    }
    std::printf("  good point n=300: h=%.1f (expected 0.0)\n", findGramOffset(300, 0.01));

    // proveSpan: true needs both endpoints good, g_m > 168*pi, and the Turing
    // bound < 2. At m=300 the bound with no bad points is 2.0204 for k=2, 1.68 for k=3.
    std::printf("\nproveSpan:\n");
    std::printf("  m=300 k=3: %d (expected 1)\n", proveSpan(300, 3, 0.01));
    std::printf("  m=300 k=2: %d (expected 0, bound 2.0204 >= 2)\n", proveSpan(300, 2, 0.01));
    std::printf("  m=365 k=3: %d (expected 1, bad point n=367 inside)\n", proveSpan(365, 3, 0.01));
    std::printf("  m=365 k=2: %d (expected 0, endpoint n=367 violates Gram's Law)\n", proveSpan(365, 2, 0.01));
    std::printf("  m=125 k=2: %d (expected 0, g_125 below 168*pi)\n", proveSpan(125, 2, 0.01));

    // findProvableSpan: smallest k in [1, kMax] that proveSpan accepts, else -1.
    std::printf("\nfindProvableSpan:\n");
    std::printf("  m=300 kMax=20: %d (expected 3)\n", findProvableSpan(300, 20, 0.01));
    std::printf("  m=365 kMax=20: %d (expected 3)\n", findProvableSpan(365, 20, 0.01));
    std::printf("  m=366 kMax=20: %d (expected 3, k=1 endpoint n=367 is bad)\n", findProvableSpan(366, 20, 0.01));
    std::printf("  m=125 kMax=20: %d (expected -1, below 168*pi)\n", findProvableSpan(125, 20, 0.01));
    std::printf("  m=300 kMax=3: %d (expected 3, kMax itself is tried)\n", findProvableSpan(300, 3, 0.01));
    std::printf("  m=300 kMax=2: %d (expected -1)\n", findProvableSpan(300, 2, 0.01));

    // proveCompleteSpan: both endpoints proven via findProvableSpan, then the
    // scan must find exactly k sign changes in between.
    std::printf("\nproveCompleteSpan:\n");
    std::printf("  m=300 k=4 kMax=20: %d (expected 1)\n", proveCompleteSpan(300, 4, 20, 0.01));
    std::printf("  m=365 k=2 kMax=20: %d (expected 0, endpoint n=367 is bad)\n", proveCompleteSpan(365, 2, 20, 0.01));
    std::printf("  m=125 k=2 kMax=20: %d (expected 0, below 168*pi)\n", proveCompleteSpan(125, 2, 20, 0.01));
    std::printf("  m=300 k=30 kMax=40: %d (expected 1, 30 zeros proven complete)\n", proveCompleteSpan(300, 30, 40, 0.01));
    std::printf("  m=300 k=0 kMax=20: %d (expected 0, k < 1)\n", proveCompleteSpan(300, 0, 20, 0.01));
    std::printf("  m=300 k=4 kMax=2: %d (expected 0, kMax too small)\n", proveCompleteSpan(300, 4, 2, 0.01));

    // Index helpers used by main: first n with g_n > 168*pi and Gram's Law
    // satisfied (n=288 is just below the threshold and bad), and the next
    // good Gram point at or after a given index (n=367 is bad).
    std::printf("\nindex helpers:\n");
    std::printf("  firstCertifiableIndex: %d (expected 289)\n", firstCertifiableIndex());
    std::printf("  nextGoodIndex(367): %d (expected 368)\n", nextGoodIndex(367));
    std::printf("  nextGoodIndex(300): %d (expected 300)\n", nextGoodIndex(300));

    // findSignChanges must not sample past tMax: a bracket straddling the end of a
    // range would be counted again by the next block's scan.
    std::printf("\nfindSignChanges end clamp:\n");
    auto clamped = findSignChanges(10.0, 55.0, 0.1);
    std::printf("  last bracket upper end: %.4f (expected <= 55.0000)\n", clamped.back().second);

    // scanBlock: a block holds nEnd - nStart zeros, and a pair closer than the step gives
    // no sign change. Block n=420889..420891 holds a pair 0.0057 apart (t ~ 273193.66),
    // missed at step 0.01 and found after one refinement. A failed call must leave `out` unchanged.
    std::printf("\nscanBlock:\n");
    std::vector<std::pair<double, double>> scanned;
    std::printf("  n=420889..420891 maxLevels=0: %d (expected 0), out size %zu (expected 0)\n",
                scanBlock(420889, 420891, 0.01, 0, scanned), scanned.size());
    std::printf("  n=420889..420891 maxLevels=4: %d (expected 1), out size %zu (expected 2)\n",
                scanBlock(420889, 420891, 0.01, 4, scanned), scanned.size());
    std::vector<std::pair<double, double>> ordinary;
    std::printf("  n=300..301 maxLevels=4: %d (expected 1), out size %zu (expected 1)\n",
                scanBlock(300, 301, 0.01, 4, ordinary), ordinary.size());

    return 0;
}
