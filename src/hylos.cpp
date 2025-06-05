#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include "tsp.h"
#include "hylos_tsp.h"
#include <cstring>

using namespace std;

int main(int argc, char* argv[]) {
    bool debug = false;
    std::string tspFile, tourFile;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--debug") == 0) {
            debug = true;
        } else if (tspFile.empty()) {
            tspFile = argv[i];
        } else if (tourFile.empty()) {
            tourFile = argv[i];
        }
    }
    
    if (tspFile.empty() || tourFile.empty()) {
        std::cerr << "Usage: " << argv[0] << " [--debug] <tsp_file> <tour_file>\n";
        return 1;
    }

    // Read coordinates
    auto coords = readTSPLib(tspFile);
    
    // Get ground truth cost
    double gtCost = tourPathCost(coords, tourFile);

    // Solve TSP
    HylosTSP solver(coords, 1.0, debug);
    
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<int> result = solver.solve();
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> elapsed = end - start;
    
    // Calculate tour cost
    double cost = pathCost(result, coords);
    
    // Calculate gap
    int gapCount = computePermutationGap(result, tourFile);
    
    // Print results using the standard format
    printResults("Hylos", cost, elapsed.count(), gtCost, gapCount);
    
    return 0;
}