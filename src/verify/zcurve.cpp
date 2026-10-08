#include <cstdio>
#include <fstream>
#include <iomanip>

#include "rs.hpp"

// Writes Z(t) from the production z(t) on a fine grid, for the first plot.
// 10 to 70 covers the first 17 zeros (14.13 ... 69.55); further out the markers start to touch.
int main() {
    // Written by Claude

    std::ofstream out("data/z_curve.csv");
    if (!out) {
        std::fprintf(stderr, "failed to open data/z_curve.csv for writing\n");
        return 1;
    }

    out << "t,z\n";
    out << std::setprecision(12);

    const double tMin = 10.0;
    const double tMax = 70.0;
    const double step = 0.01;

    long steps = static_cast<long>((tMax - tMin) / step + 0.5);
    for (long i = 0; i <= steps; ++i) {
        double t = tMin + i * step;
        out << t << "," << z(t) << "\n";
    }

    std::printf("wrote data/z_curve.csv (%ld points, t = %.0f to %.0f)\n", steps + 1, tMin, tMax);
    return 0;
}
