CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -Isrc

# All of the algorithm lives in headers, so any header change rebuilds both programs.
HEADERS = $(wildcard src/*.hpp)

# Number of zeros for `make run`. Override on the command line:
#   make run N=100
N = 20

.PHONY: all run test counts errors timing results clean

all: build/main build/tests build/counts build/errors

build/main: src/main.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) src/main.cpp -o $@

build/tests: src/verify/tests.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) src/verify/tests.cpp -o $@

build/counts: src/verify/counts.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) src/verify/counts.cpp -o $@

build/errors: src/verify/errors.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) src/verify/errors.cpp -o $@

# Find the first N zeros and write them to data/zeros.csv (overwrites it).
run: build/main
	@mkdir -p data
	./build/main $(N)

# Run the sanity-check suite. It prints results next to their expected values;
# it does not assert, so read the output.
test: build/tests
	./build/tests

# Compare zero counts with Riemann-von Mangoldt. Generates its own 200,000-zero CSV (about 30 s)
# in build/verify/, so data/zeros.csv is never touched.
counts: build/main build/counts
	@mkdir -p build/verify/data
	@cd build/verify && ../main 200000 > /dev/null && ../counts

# Error of Z(t) against the number of correction terms, against mpmath (needs `poetry install`).
errors: build/errors
	@mkdir -p data
	./build/errors
	poetry run python plot/term_errors.py

# Wall-clock time of `./build/main N` for N = 10^3 .. 10^P, run in build/verify/ so data/zeros.csv
# is never touched. Override: make timing P=6 (about 4 minutes more).
P = 5
timing: build/main
	./src/verify/timing.sh $(P)

# Every headline number in one command (about 5 minutes, mostly the 10^6-zero timing run).
# Needs `poetry install` for the mpmath comparison. Only data/term_errors.csv is overwritten:
# the counts and timing runs happen in build/verify/, so data/zeros.csv is left alone.
results: all
	@echo "== First 10 zeros: error against published values"
	@./build/tests | sed -n '/refined zeros vs known values/,/^$$/p'
	@echo "== Zero counts against Riemann-von Mangoldt (200,000 zeros)"
	@$(MAKE) --no-print-directory counts
	@echo
	@echo "== Error against number of correction terms (reference: mpmath)"
	@$(MAKE) --no-print-directory errors
	@echo
	@echo "== Timing, 10^3 to 10^6 zeros (a failed check in main stops the run)"
	@$(MAKE) --no-print-directory timing P=6

clean:
	rm -rf build
