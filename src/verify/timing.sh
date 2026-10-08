#!/bin/sh
# Wall-clock time of ./build/main N for N = 10^3 .. 10^P (default P = 5; 10^6 takes about 4 minutes).
# If time grows like T^(3/2) (about sqrt(t) terms per evaluation, about T / step evaluations),
# seconds / T^1.5 should be roughly constant.
# Runs in build/verify/, so data/zeros.csv is never touched.

P=${1:-5}

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
MAIN="$ROOT/build/main"
mkdir -p "$ROOT/build/verify/data"
cd "$ROOT/build/verify" || exit 1

printf '%9s %12s %9s %16s\n' zeros "last zero" seconds "seconds / T^1.5"

power=3
while [ "$power" -le "$P" ]; do
    n=$(awk -v p="$power" 'BEGIN { printf "%d", 10 ^ p }')

    # /usr/bin/time -p writes "real <seconds>" to the -o file. main exits non-zero when a check
    # fails (unmatched block, unproven span), so stop rather than report a time for a failed run.
    stamp=$(mktemp)
    if ! /usr/bin/time -p -o "$stamp" "$MAIN" "$n" > /dev/null; then
        echo "FAILED: $MAIN $n exited non-zero" >&2
        rm -f "$stamp"
        exit 1
    fi
    seconds=$(awk '/^real/ { print $2 }' "$stamp")
    rm -f "$stamp"
    last=$(tail -n 1 data/zeros.csv | cut -d, -f2)

    awk -v n="$n" -v t="$last" -v s="$seconds" \
        'BEGIN { printf "%9d %12.3f %9.2f %16.3e\n", n, t, s, s / (t ^ 1.5) }'

    power=$((power + 1))
done
