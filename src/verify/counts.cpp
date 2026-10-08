#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "rs.hpp"

// Reads the t column of the CSV written by main ("n,t", with a header row).
inline std::vector<double> readZeros(const std::string& path) {

    std::vector<double> zeros;
    std::ifstream in(path);

    if (!in) {
        std::fprintf(stderr, "failed to open %s for reading\n", path.c_str());
        return zeros;
    }

    std::string line;
    std::getline(in, line); // header

    while (std::getline(in, line)) {
        size_t comma = line.find(',');
        if (comma == std::string::npos) {
            continue;
        }
        zeros.push_back(std::stod(line.substr(comma + 1)));
    }

    return zeros;
}

// Number of zeros with t <= T. The CSV is in increasing order, so stop at the first larger one.
inline long countUpTo(const std::vector<double>& zeros, double T) {

    long count = 0;

    for (double t : zeros) {
        if (t > T) {
            break;
        }
        ++count;
    }

    return count;
}

int main() {

    auto zeros = readZeros("data/zeros.csv");

    double heights[] = {100.0, 1000.0, 10000.0, 100000.0};

    std::printf("%10s %10s %14s %10s %10s\n", "T", "counted", "smooth", "diff", "S(T)");

    for (double T : heights) {
        // The file must reach past T, otherwise the count is truncated.
        if (zeros.empty() || zeros.back() < T) {
            std::printf("%10.0f  skipped: data/zeros.csv does not reach this height\n", T);
            continue;
        }

        long counted = countUpTo(zeros, T);

        // Smooth part of Riemann-von Mangoldt: (T/2pi) ln(T/2pi) - T/2pi + 7/8
        double x = T / (2.0 * pi);
        double smooth = x * std::log(x) - x + 7.0 / 8.0;

        // N(T) = theta(T)/pi + 1 + S(T), so S(T) = counted - theta(T)/pi - 1
        double s = counted - theta(T) / pi - 1.0;

        std::printf("%10.0f %10ld %14.4f %10.4f %10.4f\n", T, counted, smooth, counted - smooth, s);
    }

    return 0;
}
