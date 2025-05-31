#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include "tsp.h"
#include "hylos_tsp.h"

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

    HylosTSP solver(coords, 1.2);  // alpha = 1.2
    auto start = chrono::high_resolution_clock::now();
    vector<int> path = solver.solve();
    auto end = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration_cast<chrono::duration<double>>(end - start).count();

    if (path.empty()) {
        cerr << "Hylos Failed." << endl;
        return 1;
    }

    // Save tour
    ofstream out("result/hylos_" + datasetName + ".tour");
    out << "NAME : hylos_" << datasetName << "\nTYPE : TOUR\nDIMENSION : " << coords.size() << "\nTOUR_SECTION\n";
    for (int idx : path)
        out << (idx + 1) << "\n";
    out << "-1\n";
    out.close();

    // Save clusters
    ofstream clustersOut("result/hylos_" + datasetName + "_clusters.txt");
    const auto& clusters = solver.get_clusters();
    for (size_t i = 0; i < clusters.size(); ++i) {
        clustersOut << "Cluster " << i + 1 << ":\n";
        // 센트로이드 좌표 출력
        Point centroid = solver.get_cluster_centroid(clusters[i]);
        clustersOut << "Centroid: (" << centroid.x << ", " << centroid.y << ")\n";
        // 도시 목록 출력
        clustersOut << "Cities: ";
        for (int city : clusters[i].cities) {
            clustersOut << city + 1 << " ";
        }
        clustersOut << "\n\n";
    }
    clustersOut.close();

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

    cout << "Hylos Tour cost: " << cost << endl;
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