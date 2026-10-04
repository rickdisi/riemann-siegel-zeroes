#include <cstdio>
#include <string>
#include <vector>

#include "rs.hpp"
#include "csv_export.hpp"


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: %s <targetCount> <step>\n", argv[0]);
        return 1;
    }

    // Fixed at 10.0: safely inside theta(t)'s valid range, and comfortably
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

    std::printf("found %zu zeros (target %d, scanned up to t=%.2f), written to data/zeros.csv\n",
                zeros.size(), targetCount, tMax);

    return 0;
}
