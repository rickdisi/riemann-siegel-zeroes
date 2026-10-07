CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# All of the algorithm lives in headers, so any header change rebuilds both programs.
HEADERS = $(wildcard src/*.hpp)

# Number of zeros for `make run`. Override on the command line:
#   make run N=100
N = 20

.PHONY: all run test clean

all: build/main build/tests

build/main: src/main.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) src/main.cpp -o $@

build/tests: src/tests.cpp $(HEADERS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) src/tests.cpp -o $@

# Find the first N zeros and write them to data/zeros.csv (overwrites it).
run: build/main
	@mkdir -p data
	./build/main $(N)

# Run the sanity-check suite. It prints results next to their expected values;
# it does not assert, so read the output.
test: build/tests
	./build/tests

clean:
	rm -rf build
