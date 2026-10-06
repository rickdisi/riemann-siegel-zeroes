#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>

#include "rs.hpp"
#include "turing.hpp"
#include "csv_export.hpp"


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: %s <targetCount> <step>\n", argv[0]);
        return 1;
    }

    // Fixed at 10.0: safely inside theta(t)'s valid range and comfortably
    // below the first real zero (14.134725...), so the scan always starts
    // before it and can bracket it correctly.
    double tMin = 10.0;
    int targetCount = std::stoi(argv[1]);
    double step = std::stod(argv[2]);

    double tMax = estimateTMax(targetCount);

    auto brackets = findSignChanges(tMin, tMax, step);

    std::vector<double> zeros;
    for (auto& bracket : brackets) {
        zeros.push_back(bisect(bracket.first, bracket.second, 50));
    }

    // estimateTMax() deliberately overshoots (safety margin), so trim back
    // down to exactly the number of zeros actually requested.
    if (zeros.size() > static_cast<size_t>(targetCount)) {
        zeros.resize(targetCount);
    }

    writeZerosCSV("data/zeros.csv", zeros);

    int mCert = firstCertifiableIndex();
    int nEnd = nextGoodIndex(targetCount - 1);
    int nEmpiricalEnd = std::min(mCert, nEnd); // smallest of the two

    auto failedBlocks = verifyAllGramBlocks(0, nEmpiricalEnd, step);

    std::printf("empirical Gram-block check on n=0.. %d: %zu blocks failed\n", nEmpiricalEnd, failedBlocks.size());

    if (nEnd > mCert) {
        int k = nEnd - mCert;
        int kMax = 40; // Chosen to be 40.
        bool proven = proveCompleteSpan(mCert, k, kMax, step);
        
        std::printf("Turing certification on n=%d...%d: %s\n", mCert, nEnd, proven ? "proven" : "NOT proven");
    } 
    else {
        std::printf("All requested zeroes are below Turing's minimum threshold, none certified\n");
    }

    std::printf("found %zu zeros (target %d, scanned up to t=%.2f), written to data/zeros.csv\n",
                zeros.size(), targetCount, tMax);

    return 0;
}
