CXX      = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Compile and run the algorithm.
# Usage: make run FILE=src/rs.cpp
run:
	$(CXX) $(CXXFLAGS) $(FILE) -o /tmp/rsrun && /tmp/rsrun
