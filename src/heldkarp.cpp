#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include "tsp.h"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <tour_file>" << endl;
        return 1;
    }

    string tspFile = argv[1];
    string gtTourFile = argv[2];
    string datasetName = tspFile.substr(tspFile.find_last_of("/") + 1);
    datasetName = datasetName.substr(0, datasetName.find_last_of("."));

    vector<Point> coords = readTSPLib(tspFile);
    if (coords.empty()) {
        cerr << "Failed to load coordinates." << endl;
        return 1;
    }

    if (coords.size() > 22) {
        cerr << "Warning: Held-Karp is impractical for n > 22" << endl;
        cerr << "Current size: " << coords.size() << endl;
        return 1;
    }

    auto start = chrono::high_resolution_clock::now();
    vector<int> path = heldkarpPath(coords);
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration_cast<chrono::duration<double>>(end - start).count();

    if (path.empty()) {
        cerr << "Held-Karp Failed." << endl;
        return 1;
    }

    // Save tour
    ofstream out("result/heldkarp_" + datasetName + ".tour");
    if (!out.is_open()) {
        cerr << "Failed to create output file" << endl;
        return 1;
    }
    out << "NAME : heldkarp_" << datasetName << "\nTYPE : TOUR\nDIMENSION : " << coords.size() << "\nTOUR_SECTION\n";
    for (int idx : path)
        out << (idx + 1) << "\n";
    out << "-1\n";
    out.close();

    // Compute cost
    double cost = pathCost(path, coords);
    if (cost < 0) {
        cerr << "Invalid path cost" << endl;
        return 1;
    }

    // Compare with ground truth tour
    double gtCost = tourPathCost(coords, gtTourFile);

    // Compute gap
    int gapCount = computePermutationGap(path, gtTourFile);

    printResults("Held-Karp", cost, elapsed, gtCost, gapCount);

    return 0;
}