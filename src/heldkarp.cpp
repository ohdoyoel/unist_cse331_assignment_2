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
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <dataset.tsp>" << endl;
        return 1;
    }

    string tspFile = argv[1];
    string datasetName = tspFile.substr(tspFile.find_last_of("/") + 1);
    datasetName = datasetName.substr(0, datasetName.find_last_of("."));

    vector<Point> coords = readTSPLib(tspFile);
    if (coords.empty()) {
        cerr << "Failed to load coordinates." << endl;
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
    out << "NAME : heldkarp_" << datasetName << "\nTYPE : TOUR\nDIMENSION : " << coords.size() << "\nTOUR_SECTION\n";
    for (int idx : path)
        out << (idx + 1) << "\n";
    out << "-1\n";
    out.close();

    // Compute cost
    double cost = pathCost(path, coords);

    // Compare with ground truth tour
    string gtTourFile = "data/" + datasetName + ".tour";
    double gtCost = tourPathCost(coords, gtTourFile);

    // Compute gap
    int gapCount = 0;
    if (gtCost > 0) {
        gapCount = computePermutationGap(path, gtTourFile);
    }

    cout << "Held-Karp Tour cost: " << cost << endl;
    cout << "Elapsed time: " << elapsed << " sec" << endl;
    if (gtCost > 0)
        cout << "Optimal Tour cost: " << gtCost 
        << "\nApproximation Quality: " << (gtCost / cost * 100.0) << "%"
        << "\nGap: " << gapCount << " / " << coords.size()
        << endl;
    else
        cout << "Ground-truth tour not found or invalid." << endl;

    return 0;
}