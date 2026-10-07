#include <cstdio>
#include <string>
#include <vector>

#include "rs.hpp"
#include "turing.hpp"
#include "csv_export.hpp"


int main(int argc, char* argv[]) {
    
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s <targetCount>\n", argv[0]);
        return 1;
    }

    // Longest span findProvableSpan will try when proving a zero count at one Gram point.
    const int kMax = 40;
    // Scan step in t. A pair of zeros closer than this is missed at the first pass, then caught by the per-block rescan below.
    const double baseStep = 0.01;
    const int maxLevels = 4; // Each rescan divides the step by 10: hence 0.01 down to 1e-6 at maxLevels = 4.

    // Fixed at 10.0: safely inside theta(t)'s valid range and below the first zero (14.134725...).
    double tMin = 10.0;

    int targetCount = std::stoi(argv[1]);
    int mCert = firstCertifiableIndex(); // First Gram point above 168*pi: Turing's bound is valid from here on.
    int nEnd = nextGoodIndex(targetCount - 1); // Last Gram point of the run.

    std::vector<std::pair<double, double>> allBrackets;

    // Zero 1 lies below g_0, so it isn't inside any Gram block
    auto firstBracket = findSignChanges(tMin, gram(0), baseStep);
    if (firstBracket.size() != 1) {
        std::fprintf(stderr, "FAILED: Could not find zero 1 below g_0.\n");
        return 1;
    }
    allBrackets.insert(allBrackets.end(), firstBracket.begin(), firstBracket.end());

    // Brackets found at or above mCert, to compare against the proven zero count.
    long certifiedCount = 0;
    auto blocks = findGramBlocks(0, nEnd); // Consecutive good Gram points bound each block, each holding a known number of zeros.

    for (auto& block : blocks) {
        size_t before = allBrackets.size();

        if (!scanBlock(block.first, block.second, baseStep, maxLevels, allBrackets)) {
            std::fprintf(stderr, "FAILED: Gram block n=%d..%d did not match its expected zero count at any step\n", block.first, block.second);
            return 1;
        }
        if (block.first >= mCert) {
            certifiedCount += static_cast<long>(allBrackets.size() - before);
        }
    }

    // Refine every bracket to a zero location.
    std::vector<double> zeros;
    for (auto& bracket : allBrackets) {
        zeros.push_back(bisect(bracket.first, bracket.second, 50));
    }

    // The span ends at a Gram point, so it may hold a few zeros past the target.
    if (zeros.size() > static_cast<size_t>(targetCount)) {
        zeros.resize(targetCount);
    }

    // Above Turing's threshold the zero count is proven, not just checked:
    // N is proven at both ends of the span, and the scan must have found exactly that many zeros in between
    if (nEnd > mCert) {
        int k = nEnd - mCert;
        bool endsProven = findProvableSpan(mCert, kMax, baseStep) != -1 && findProvableSpan(nEnd, kMax, baseStep) != -1;
        bool countOk = certifiedCount == k;

        if (!endsProven) {
            std::fprintf(stderr, "FAILED: could not prove the zero counts at n=%d or n=%d\n", mCert, nEnd);
            return 1;
        }
        if (!countOk) {
            std::fprintf(stderr, "FAILED: found %ld zeros in n=%d..%d, expected %d\n", certifiedCount, mCert, nEnd, k);
            return 1;
        }
        std::printf("Turing certification on n=%d..%d: proven (%d zeros)\n", mCert, nEnd, k);
    } else {
        std::printf("No requested zero is above Turing's threshold; Gram blocks matched their expected counts at every level\n");
    }

    // Written last, so a failed check above never leaves a CSV behind.
    writeZerosCSV("data/zeros.csv", zeros);

    std::printf("Found %zu zeros up to Gram point n=%d, written to data/zeros.csv\n", zeros.size(), nEnd);

    return 0;
}