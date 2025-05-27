#include "tsp.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stack>
#include <set>
#include <limits>
#include <algorithm>

using namespace std;

vector<Point> readTSPLib(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Cannot open file: " << filename << endl;
        return {};
    }

    string line;
    vector<Point> coords;
    bool start = false;
    while (getline(file, line)) {
        if (line.find("NODE_COORD_SECTION") != string::npos) {
            start = true;
            continue;
        }
        if (line.find("EOF") != string::npos) break;
        if (start) {
            stringstream ss(line);
            int idx;
            double x, y;
            ss >> idx >> x >> y;
            if (ss) coords.push_back({x, y});
        }
    }

    return coords;
}

double euclideanDistance(const Point& a, const Point& b) {
    return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

// vector<vector<double>> computeDistanceMatrix(const vector<Point>& coords) {
//     int n = coords.size();
//     vector<vector<double>> dist(n, vector<double>(n));
//     for (int i = 0; i < n; ++i)
//         for (int j = 0; j < n; ++j)
//             dist[i][j] = euclideanDistance(coords[i], coords[j]);
//     return dist;
// }

// Prim's Algorithm for MST
vector<pair<int, int>> minimumSpanningTree(const vector<Point>& coords) {
    int n = coords.size();
    vector<bool> inTree(n, false);
    vector<double> minEdge(n, numeric_limits<double>::max());
    vector<int> parent(n, -1);
    vector<pair<int, int>> mst;

    minEdge[0] = 0;

    for (int i = 0; i < n; ++i) {
        double minCost = numeric_limits<double>::max();
        int u = -1;

        for (int j = 0; j < n; ++j) {
            if (!inTree[j] && minEdge[j] < minCost) {
                minCost = minEdge[j];
                u = j;
            }
        }

        inTree[u] = true;
        if (parent[u] != -1) mst.emplace_back(u, parent[u]);

        for (int v = 0; v < n; ++v) {
            double d = euclideanDistance(coords[u], coords[v]);
            if (!inTree[v] && d < minEdge[v]) {
                minEdge[v] = d;
                parent[v] = u;
            }
        }
    }

    return mst;
}

// Naive matching: greedy pair closest odd-degree nodes (for large data use)
vector<pair<int, int>> greedyPerfectMatching(const vector<int>& odd, const vector<Point>& coords) {
    vector<pair<int, int>> matching;
    vector<bool> used(odd.size(), false);
    for (int i = 0; i < odd.size(); ++i) {
        if (used[i]) continue;
        int best = -1;
        double bestDist = numeric_limits<double>::max();
        for (int j = i + 1; j < odd.size(); ++j) {
            double d = euclideanDistance(coords[odd[i]], coords[odd[j]]);
            if (!used[j] && d < bestDist) {
                best = j;
                bestDist = d;
            }
        }
        if (best != -1) {
            matching.emplace_back(odd[i], odd[best]);
            used[i] = used[best] = true;
        }
    }
    return matching;
}

// Euler tour using Hierholzer's Algorithm
vector<int> eulerianTour(multimap<int, int>& graph, int start) {
    vector<int> tour;
    stack<int> s;
    s.push(start);

    while (!s.empty()) {
        int u = s.top();
        auto range = graph.equal_range(u);
        if (range.first == range.second) {
            tour.push_back(u);
            s.pop();
        } else {
            int v = range.first->second;
            graph.erase(range.first);
            auto it = graph.find(v);
            while (it != graph.end() && it->second != u) ++it;
            if (it != graph.end()) graph.erase(it);
            s.push(v);
        }
    }
    return tour;
}

vector<int> makeHamiltonian(const vector<int>& eulerTour) {
    set<int> visited;
    vector<int> path;
    for (int node : eulerTour) {
        if (!visited.count(node)) {
            visited.insert(node);
            path.push_back(node);
        }
    }
    return path;
}

double pathCost(const vector<int>& path, const vector<Point>& coords) {
    double cost = 0;
    for (size_t i = 0; i + 1 < path.size(); ++i) {
        cost += euclideanDistance(coords[path[i]], coords[path[i + 1]]);
    }
    return cost;
}

vector<int> parseTourFile(const string& tourFile) {
    ifstream in(tourFile);
    if (!in.is_open()) {
        cerr << "Cannot open ground-truth tour file: " << tourFile << endl;
        return {};
    }

    string line;
    while (getline(in, line)) {
        if (line.find("TOUR_SECTION") != string::npos)
            break;
    }

    vector<int> tour;
    while (getline(in, line)) {
        istringstream iss(line);
        int idx;
        while (iss >> idx) {
            if (idx == -1) return tour;
            tour.push_back(idx - 1);  // 0-based
        }
    }
    in.close();

    return tour;
}

double tourPathCost(const vector<Point>& coords, const string& tourFile) {
    vector<int> gtPath = parseTourFile(tourFile);

    if (gtPath.size() != coords.size()) {
        cerr << "tourPathCost() : Mismatch in tour length: expected " << coords.size() << ", got " << gtPath.size() << endl;
        return -1;
    }

    return pathCost(gtPath, coords);
}

int computePermutationGap(const vector<int>& result, const string& tourFile) {
    vector<int> optimal = parseTourFile(tourFile);

    if (optimal.size() != result.size()) {
        cerr << "computePermutationGap() : Mismatch in tour length: expected " << optimal.size() << ", got " << result.size() << endl;
        return -1;
    }
    
    int n = result.size();
    int minGap = n;

    for (int offset = 0; offset < n; ++offset) {
        int mismatch = 0;
        for (int i = 0; i < n; ++i) {
            if (result[i] != optimal[(i + offset) % n])
                mismatch++;
        }
        minGap = min(minGap, mismatch);
    }

    vector<int> reversed = result;
    reverse(reversed.begin(), reversed.end());
    for (int offset = 0; offset < n; ++offset) {
        int mismatch = 0;
        for (int i = 0; i < n; ++i) {
            if (reversed[i] != optimal[(i + offset) % n])
                mismatch++;
        }
        minGap = min(minGap, mismatch);
    }

    return minGap;
}


// Christofides Algorithm for TSP
vector<int> christofidesPath(const vector<Point>& coords) {
    auto mst = minimumSpanningTree(coords);

    map<int, int> degree;
    for (auto& e : mst) {
        degree[e.first]++;
        degree[e.second]++;
    }

    vector<int> odd;
    for (auto& [v, d] : degree)
        if (d % 2 == 1) odd.push_back(v);

    auto matching = greedyPerfectMatching(odd, coords);

    multimap<int, int> multigraph;
    for (auto& [u, v] : mst) {
        multigraph.insert({u, v});
        multigraph.insert({v, u});
    }
    for (auto& [u, v] : matching) {
        multigraph.insert({u, v});
        multigraph.insert({v, u});
    }

    auto euler = eulerianTour(multigraph, 0);
    return makeHamiltonian(euler);
}

// Held-Karp Algorithm for TSP
vector<int> heldkarpPath(const vector<Point>& coords) {
    const int n = coords.size();
    // if (n > 30) {
    //     cerr << "Held-Karp is too slow for n > 30." << endl;
    //     return {};
    // }

    const double INF = numeric_limits<double>::max() / 2;
    vector<vector<double>> dist(n, vector<double>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            dist[i][j] = euclideanDistance(coords[i], coords[j]);

    vector<vector<double>> dp(1 << n, vector<double>(n, INF));
    vector<vector<int>> parent(1 << n, vector<int>(n, -1));
    dp[1][0] = 0;

    for (int mask = 1; mask < (1 << n); ++mask) {
        for (int u = 0; u < n; ++u) {
            if (!(mask & (1 << u))) continue;
            for (int v = 0; v < n; ++v) {
                if (mask & (1 << v)) continue;
                int next = mask | (1 << v);
                double newCost = dp[mask][u] + dist[u][v];
                if (newCost < dp[next][v]) {
                    dp[next][v] = newCost;
                    parent[next][v] = u;
                }
            }
        }
    }

    // Reconstruct path
    int last = -1;
    double minCost = INF;
    int full = (1 << n) - 1;
    for (int i = 1; i < n; ++i) {
        double cost = dp[full][i] + dist[i][0];
        if (cost < minCost) {
            minCost = cost;
            last = i;
        }
    }

    vector<int> path;
    int mask = full;
    while (last != -1) {
        path.push_back(last);
        int temp = parent[mask][last];
        mask ^= (1 << last);
        last = temp;
    }

    reverse(path.begin(), path.end());
    return path;
}
