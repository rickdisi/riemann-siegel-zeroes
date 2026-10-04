#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

// Writes the refined zero locations to a CSV file at `path`, one per line.
inline void writeZerosCSV(const std::string& path, const std::vector<double>& zeros) {

    std::ofstream out(path);

    if (!out) {
        std::fprintf(stderr, "failed to open %s for writing\n", path.c_str());
        return;
    }
    out << "n,t\n";

    out << std::setprecision(12); 
    for (size_t n = 1; n <= zeros.size(); ++n) {
        out << n << "," << zeros[n - 1] << "\n";
    }
}