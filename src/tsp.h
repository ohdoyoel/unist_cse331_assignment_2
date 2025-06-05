#ifndef TSP_H
#define TSP_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>

struct Point {
    double x, y;
};

// Common output function for all algorithms
inline void printResults(const std::string& algorithm_name, double cost, double elapsed, double gtCost, int gapCount) {
    std::cout << "Tour cost: " << cost << "\n"
              << "Elapsed time: " << elapsed << " sec\n";
    if (gtCost > 0) {
        std::cout << "Optimal Tour cost: " << gtCost << "\n"
                  << "Approximation Ratio: " << (cost / gtCost) << "\n"
                  << "Gap: " << gapCount;
    }
}

double euclideanDistance(const Point& a, const Point& b);
std::vector<Point> readTSPLib(const std::string& filename);
std::vector<std::vector<double>> computeDistanceMatrix(const std::vector<Point>& coords);
std::vector<std::pair<int, int>> minimumSpanningTree(const std::vector<Point>& coords);
std::vector<std::pair<int, int>> greedyPerfectMatching(const std::vector<int>& odd, const std::vector<Point>& coords);
std::vector<int> eulerianTour(std::multimap<int, int>& graph, int start);
std::vector<int> makeHamiltonian(const std::vector<int>& eulerTour);
double pathCost(const std::vector<int>& path, const std::vector<Point>& coords);
double tourPathCost(const std::vector<Point>& coords, const std::string& tourFile);
int computePermutationGap(const std::vector<int>& result, const std::string& tourFile);
std::vector<int> christofidesPath(const std::vector<Point>& coords);
std::vector<int> christofidesPath_edmonds(const std::vector<Point>& coords);
std::vector<int> christofidesPath_gabow(const std::vector<Point>& coords);
std::vector<int> heldkarpPath(const std::vector<Point>& coords);
std::vector<int> mst2approx(const std::vector<Point>& coords);

#endif