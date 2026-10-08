"""
Error of Z(t) against the number of Riemann-Siegel correction terms.

Reads data/term_errors.csv (written by src/verify/errors.cpp), evaluates the reference
Z(t) with mpmath at 25 digits, and prints the maximum error over p at each
height for C0 only, C0..C1 and C0..C2, plus the fitted log-log slope.

The first omitted term is C_{K+1}(p) (2 pi / t)^(1/4 + (K+1)/2), so the error
should scale as t^(-3/4), t^(-5/4), t^(-7/4) for K = 0, 1, 2. mpmath uses its own
Riemann-Siegel routine at large t, so it is independent of our truncation,
coefficients and double precision, not of the method itself.
"""

import csv
from collections import defaultdict

import mpmath
import numpy as np

mpmath.mp.dps = 25

PREDICTED = {1: -0.75, 2: -1.25, 3: -1.75}

# Below this the error is double-precision rounding of the phase (it grows with t), not
# truncation, so those heights are left out of the slope fit. Chosen after looking at
# where the C0..C2 error stops falling (~3e-10 near t = 2e4); it is not derived.
FLOOR = 5e-9

# max |error| over p, per height N, per number of terms
worst = defaultdict(lambda: {1: 0.0, 2: 0.0, 3: 0.0})
height = {}

with open("data/term_errors.csv") as f:
    for row in csv.DictReader(f):
        n = int(row["N"])
        t = float(row["t"])
        ref = float(mpmath.siegelz(mpmath.mpf(t)))
        height[n] = min(height.get(n, t), t)
        for k, col in ((1, "z1"), (2, "z2"), (3, "z3")):
            worst[n][k] = max(worst[n][k], abs(float(row[col]) - ref))

ns = sorted(worst)
print(f"{'N':>5} {'t (p~0)':>12} {'C0 only':>12} {'C0..C1':>12} {'C0..C2':>12}")
for n in ns:
    print(f"{n:>5} {height[n]:>12.1f} {worst[n][1]:>12.3e} {worst[n][2]:>12.3e} {worst[n][3]:>12.3e}")

print()
for k, label in ((1, "C0 only"), (2, "C0..C1"), (3, "C0..C2")):
    used = [n for n in ns if worst[n][k] > FLOOR]
    x = np.log([height[n] for n in used])
    y = np.log([worst[n][k] for n in used])
    slope = np.polyfit(x, y, 1)[0]
    print(f"{label:>8}: fitted slope {slope:+.2f} over {len(used)} heights with error > {FLOOR:.0e}, predicted {PREDICTED[k]:+.2f}")
