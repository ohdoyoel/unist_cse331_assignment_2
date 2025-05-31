#ifndef TSP_H
#define TSP_H

#include <vector>
#include <string>
#include <map>

struct Point {
    double x, y;
};

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
std::vector<int> heldkarpPath(const std::vector<Point>& coords);
std::vector<int> mst2approx(const std::vector<Point>& coords);

#endif