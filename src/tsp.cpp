#include "tsp.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stack>
#include <set>
#include <limits>
#include <algorithm>
#include <queue>
#include <map>
#include <unordered_map>
#include <functional>

using namespace std;

double euclideanDistance(const Point& a, const Point& b) {
    return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

// Blossom algorithm data structures
struct BlossomNode {
    int parent;
    int mate;
    int label;  // 0: unlabeled, 1: outer, -1: inner
    int blossom_root;  // root of the blossom containing this vertex
    double dual;  // dual variable for primal-dual method
    
    BlossomNode() : parent(-1), mate(-1), label(0), blossom_root(-1), dual(0) {}
};

struct Edge {
    int u, v;
    double weight;
    Edge(int _u, int _v, double _w) : u(_u), v(_v), weight(_w) {}
};

// Helper functions for Blossom algorithms
int find_blossom_root(vector<BlossomNode>& nodes, int v) {
    if (v < 0 || v >= nodes.size()) return -1;  // Safety check
    if (nodes[v].blossom_root == -1) return v;
    nodes[v].blossom_root = find_blossom_root(nodes, nodes[v].blossom_root);
    return nodes[v].blossom_root;
}

vector<int> find_augmenting_path(vector<BlossomNode>& nodes, int start, int end) {
    vector<int> path;
    if (start < 0 || end < 0 || start >= nodes.size() || end >= nodes.size()) {
        return path;  // Return empty path if invalid indices
    }
    
    vector<bool> visited(nodes.size(), false);  // Prevent infinite loops
    int current = end;
    while (current != start && current != -1 && !visited[current]) {
        visited[current] = true;
        path.push_back(current);
        current = nodes[current].parent;
    }
    
    if (current == start) {
        path.push_back(start);
        reverse(path.begin(), path.end());
        return path;
    }
    
    return vector<int>();  // Return empty path if no valid path found
}

// Helper function to get sorted edges
vector<Edge> getSortedEdges(const vector<Point>& coords) {
    vector<Edge> edges;
    int n = coords.size();
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double weight = euclideanDistance(coords[i], coords[j]);
            edges.emplace_back(i, j, weight);
        }
    }
    sort(edges.begin(), edges.end(), 
         [](const Edge& a, const Edge& b) { return a.weight < b.weight; });
    return edges;
}

// Edmonds' Blossom Algorithm O(n^4)
vector<pair<int, int>> blossom_edmonds(const vector<Point>& coords) {
    if (coords.empty()) return {};
    
    int n = coords.size();
    vector<BlossomNode> nodes(n);
    vector<Edge> edges = getSortedEdges(coords);  // Use sorted edges
    
    vector<pair<int, int>> matching;
    
    // Main algorithm loop
    int max_iterations = n * n;  // Prevent infinite loops
    int iteration = 0;
    
    while (iteration++ < max_iterations) {
        // Reset labels and find an exposed vertex
        for (auto& node : nodes) {
            node.label = 0;
            node.parent = -1;
            node.blossom_root = -1;  // Reset blossom root
        }
        
        int exposed = -1;
        for (int i = 0; i < n; i++) {
            if (nodes[i].mate == -1) {
                exposed = i;
                break;
            }
        }
        
        if (exposed == -1) break;  // Perfect matching found
        
        // BFS to find augmenting path
        queue<int> q;
        nodes[exposed].label = 1;  // Mark as outer
        q.push(exposed);
        
        bool found_augmenting_path = false;
        while (!q.empty() && !found_augmenting_path) {
            int u = q.front();
            q.pop();
            
            for (const Edge& e : edges) {  // Process edges in sorted order
                // Fix the bug in vertex selection
                int v = (e.u == u) ? e.v : (e.v == u ? e.u : -1);
                if (v == -1) continue;
                
                int v_root = find_blossom_root(nodes, v);
                int u_root = find_blossom_root(nodes, u);
                
                if (v_root == -1 || u_root == -1 || v_root == u_root) continue;
                
                if (nodes[v].label == 0) {  // v is unlabeled
                    if (nodes[v].mate != -1) {
                        nodes[v].label = -1;  // inner
                        nodes[v].parent = u;
                        int mate = nodes[v].mate;
                        if (mate >= 0 && mate < n) {  // Safety check
                            nodes[mate].label = 1;  // outer
                            nodes[mate].parent = v;
                            q.push(mate);
                        }
                    } else {
                        // Found augmenting path
                        nodes[v].parent = u;
                        vector<int> path = find_augmenting_path(nodes, exposed, v);
                        
                        if (!path.empty()) {  // Only augment if valid path found
                            // Augment matching along path
                            for (size_t i = 0; i < path.size() - 1; i += 2) {
                                nodes[path[i]].mate = path[i + 1];
                                nodes[path[i + 1]].mate = path[i];
                            }
                            found_augmenting_path = true;
                            break;
                        }
                    }
                }
            }
        }
        
        if (!found_augmenting_path) break;
    }
    
    // Construct final matching
    for (int i = 0; i < n; i++) {
        if (i < nodes[i].mate) {
            matching.emplace_back(i, nodes[i].mate);
        }
    }
    
    return matching;
}

// Gabow's improved Blossom Algorithm O(n^3)
vector<pair<int, int>> blossom_gabow(const vector<Point>& coords) {
    if (coords.empty()) {
        return {};
    }
    
    int n = coords.size();
    vector<BlossomNode> nodes(n);
    vector<vector<double>> weights(n, vector<double>(n));
    
    // Initialize weight matrix
    try {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                weights[i][j] = euclideanDistance(coords[i], coords[j]);
            }
        }
    } catch (const exception& e) {
        return {};
    }
    
    // Initialize dual variables
    for (int i = 0; i < n; i++) {
        double min_weight = numeric_limits<double>::max();
        for (int j = 0; j < n; j++) {
            if (i != j) {
                min_weight = min(min_weight, weights[i][j]);
            }
        }
        nodes[i].dual = min_weight / 2;
    }
    
    vector<pair<int, int>> matching;
    int iterations = 0;
    const int MAX_ITERATIONS = n * n * n;  // O(n^3) bound
    
    while (iterations++ < MAX_ITERATIONS) {
        // Reset labels and find an exposed vertex
        for (auto& node : nodes) {
            node.label = 0;
            node.parent = -1;
            node.blossom_root = -1;
        }
        
        int exposed = -1;
        for (int i = 0; i < n; i++) {
            if (nodes[i].mate == -1) {
                exposed = i;
                break;
            }
        }
        
        if (exposed == -1) {
            break;  // Perfect matching found
        }
        
        // Priority queue for efficient edge selection
        using pqtype = tuple<double, int, int>;  // (slack, u, v)
        priority_queue<pqtype, vector<pqtype>, greater<pqtype>> pq;
        
        nodes[exposed].label = 1;
        
        // Initialize priority queue with edges from exposed vertex
        for (int v = 0; v < n; v++) {
            if (v != exposed) {
                double slack = weights[exposed][v] - nodes[exposed].dual - nodes[v].dual;
                pq.push({slack, exposed, v});
            }
        }
        
        bool augmented = false;
        set<pair<int, int>> processed_edges;  // Track processed edges to avoid cycles
        
        while (!pq.empty() && !augmented) {
            auto [slack, u, v] = pq.top();
            pq.pop();
            
            // Skip if we've already processed this edge
            if (processed_edges.count({u, v}) || processed_edges.count({v, u})) {
                continue;
            }
            processed_edges.insert({u, v});
            
            if (nodes[v].label == 0) {  // v is unlabeled
                if (nodes[v].mate == -1) {
                    // Found augmenting path
                    nodes[v].parent = u;
                    vector<int> path = find_augmenting_path(nodes, exposed, v);
                    
                    if (!path.empty()) {
                        // Augment matching along path
                        for (size_t i = 0; i < path.size() - 1; i += 2) {
                            if (path[i] >= 0 && path[i] < n && path[i+1] >= 0 && path[i+1] < n) {
                                nodes[path[i]].mate = path[i + 1];
                                nodes[path[i + 1]].mate = path[i];
                            } else {
                                return {};
                            }
                        }
                        augmented = true;
                    }
                } else {
                    // Extend alternating path
                    int mate = nodes[v].mate;
                    if (mate >= 0 && mate < n) {
                        nodes[v].label = -1;
                        nodes[v].parent = u;
                        nodes[mate].label = 1;
                        nodes[mate].parent = v;
                        
                        // Add edges from new outer vertex
                        for (int w = 0; w < n; w++) {
                            if (nodes[w].label == 0) {
                                double new_slack = weights[mate][w] - nodes[mate].dual - nodes[w].dual;
                                pq.push({new_slack, mate, w});
                            }
                        }
                    }
                }
            }
        }
        
        if (!augmented) {
            // Update dual variables
            vector<double> delta(n, numeric_limits<double>::max());
            for (int i = 0; i < n; i++) {
                if (nodes[i].label == 1) {  // outer
                    for (int j = 0; j < n; j++) {
                        if (nodes[j].label == 0) {  // unlabeled
                            delta[i] = min(delta[i], (weights[i][j] - nodes[i].dual - nodes[j].dual) / 2);
                        }
                    }
                }
            }
            
            double min_delta = numeric_limits<double>::max();
            for (double d : delta) {
                if (d > 0) min_delta = min(min_delta, d);
            }
            
            if (min_delta == numeric_limits<double>::max()) {
                break;
            }
            
            // Apply dual adjustments
            for (int i = 0; i < n; i++) {
                if (nodes[i].label == 1) nodes[i].dual += min_delta;
                else if (nodes[i].label == -1) nodes[i].dual -= min_delta;
            }
        }
    }
    
    if (iterations >= MAX_ITERATIONS) {
        return {};
    }
    
    // Construct final matching
    for (int i = 0; i < n; i++) {
        if (i < nodes[i].mate) {
            matching.emplace_back(i, nodes[i].mate);
        }
    }
    
    return matching;
}

vector<Point> readTSPLib(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
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

// Naive matching: greedy pair closest odd-degree nodes
vector<pair<int, int>> greedyPerfectMatching(const vector<int>& odd_vertices, const vector<Point>& coords) {
    if (odd_vertices.empty()) return {};
    
    vector<pair<int, int>> matching;
    vector<bool> used(odd_vertices.size(), false);

    // For each unmatched vertex, find the closest unmatched neighbor
    for (size_t i = 0; i < odd_vertices.size(); ++i) {
        if (used[i]) continue;

        int best_j = -1;
        double best_dist = numeric_limits<double>::max();

        // Find closest unmatched vertex
        for (size_t j = i + 1; j < odd_vertices.size(); ++j) {
            if (!used[j]) {
                double dist = euclideanDistance(coords[odd_vertices[i]], coords[odd_vertices[j]]);
                if (dist < best_dist) {
                    best_dist = dist;
                    best_j = j;
                }
            }
        }

        if (best_j != -1) {
            // Store the original vertex indices from odd_vertices
            matching.emplace_back(odd_vertices[i], odd_vertices[best_j]);
            used[i] = used[best_j] = true;
        }
    }

    return matching;
}

// Euler tour using Hierholzer's Algorithm
vector<int> eulerianTour(multimap<int, int>& graph, int start) {
    if (graph.empty()) {
        return {};
    }

    vector<int> tour;
    stack<int> s;
    set<int> vertices;  // Track all vertices
    
    // First, collect all vertices
    for (const auto& [u, v] : graph) {
        vertices.insert(u);
        vertices.insert(v);
    }
    
    // Verify start vertex exists in graph
    if (vertices.find(start) == vertices.end()) {
        return {};
    }
    
    // Track the number of edges to ensure we don't get stuck
    size_t total_edges = graph.size();
    size_t edges_used = 0;
    
    s.push(start);
    
    while (!s.empty()) {
        int u = s.top();
        
        auto range = graph.equal_range(u);
        if (range.first == range.second) {
            // No more edges from this vertex
            tour.push_back(u);
            s.pop();
        } else {
            // Take the next available edge
            int v = range.first->second;
            
            // Remove both directions of the edge
            graph.erase(range.first);
            edges_used++;
            
            // Find and remove the reverse edge
            auto rev_range = graph.equal_range(v);
            for (auto it = rev_range.first; it != rev_range.second; ++it) {
                if (it->second == u) {
                    graph.erase(it);
                    edges_used++;
                    break;
                }
            }
            
            s.push(v);
        }
        
        // Safety check: ensure we're not stuck
        if (edges_used > total_edges) {
            return {};
        }
    }
    
    // Verify the tour contains all vertices
    set<int> tour_vertices(tour.begin(), tour.end());
    if (tour_vertices != vertices) {
        return {};
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
    cost += euclideanDistance(coords[path[path.size() - 1]], coords[path[0]]);
    return cost;
}

vector<int> parseTourFile(const string& tourFile) {
    ifstream in(tourFile);
    if (!in.is_open()) {
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
        return -1;
    }

    return pathCost(gtPath, coords);
}

int computePermutationGap(const vector<int>& result, const string& tourFile) {
    vector<int> optimal = parseTourFile(tourFile);

    if (optimal.size() != result.size()) {
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
    if (coords.empty()) {
        return {};
    }
    
    // Get MST
    auto mst = minimumSpanningTree(coords);
    if (mst.empty()) {
        return vector<int>(coords.size());
    }

    // Find vertices with odd degree
    map<int, int> degree;
    for (auto& e : mst) {
        if (e.first >= coords.size() || e.second >= coords.size()) {
            return vector<int>(coords.size());
        }
        degree[e.first]++;
        degree[e.second]++;
    }

    vector<int> odd;
    for (auto& [v, d] : degree) {
        if (d % 2 == 1) {
            if (v >= 0 && v < coords.size()) {
                odd.push_back(v);
            } else {
                return vector<int>(coords.size());
            }
        }
    }

    // Choose matching algorithm based on problem size
    vector<pair<int, int>> matching;
    if (odd.size() > 10000) {
        matching = greedyPerfectMatching(odd, coords);
    } else {
        // Create subgraph for Gabow's algorithm
        vector<Point> odd_vertices;
        for (int idx : odd) {
            odd_vertices.push_back(coords[idx]);
        }
        
        auto gabow_matching = blossom_gabow(odd_vertices);
        
        // Convert Gabow's matching indices back to original graph indices
        for (auto [u, v] : gabow_matching) {
            if (u < odd.size() && v < odd.size()) {
                matching.emplace_back(odd[u], odd[v]);
            } else {
                return vector<int>(coords.size());
            }
        }
    }

    // Combine MST and matching edges into multigraph
    multimap<int, int> multigraph;
    try {
        for (auto& [u, v] : mst) {
            if (u >= coords.size() || v >= coords.size()) continue;
            multigraph.insert({u, v});
            multigraph.insert({v, u});
        }
        for (auto& [u, v] : matching) {
            if (u >= coords.size() || v >= coords.size()) continue;
            multigraph.insert({u, v});
            multigraph.insert({v, u});
        }
    } catch (const exception& e) {
        return vector<int>(coords.size());
    }

    if (multigraph.empty()) {
        return vector<int>(coords.size());
    }

    // Find Eulerian tour
    try {
        auto euler = eulerianTour(multigraph, 0);
        if (euler.empty()) {
            return vector<int>(coords.size());
        }
        
        // Make Hamiltonian
        auto hamilton = makeHamiltonian(euler);
        if (hamilton.empty()) {
            return vector<int>(coords.size());
        }
        
        // Verify path
        for (int v : hamilton) {
            if (v < 0 || v >= coords.size()) {
                return vector<int>(coords.size());
            }
        }
        
        return hamilton;
    } catch (const exception& e) {
        return vector<int>(coords.size());
    }
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

// MST-based 2-approximation algorithm for TSP
vector<int> mst2approx(const vector<Point>& coords) {
    if (coords.empty()) return {};
    
    int n = coords.size();
    if (n <= 1) return {0};
    
    // Get the minimum spanning tree edges
    vector<pair<int, int>> mst_edges = minimumSpanningTree(coords);
    
    // Create adjacency list representation for the MST
    vector<vector<int>> adj(n);
    for (const auto& edge : mst_edges) {
        adj[edge.first].push_back(edge.second);
        adj[edge.second].push_back(edge.first);
    }
    
    // Perform preorder traversal of the MST
    vector<bool> visited(n, false);
    vector<int> tour;
    
    std::function<void(int)> preorder = [&](int v) {
        visited[v] = true;
        tour.push_back(v);
        
        for (int u : adj[v]) {
            if (!visited[u]) {
                preorder(u);
            }
        }
    };
    
    // Start DFS from vertex 0
    preorder(0);
    
    // Add the starting vertex to complete the tour
    // tour.push_back(0);
    
    return tour;
}

// Christofides Algorithm with Edmonds' blossom matching
vector<int> christofidesPath_edmonds(const vector<Point>& coords) {
    if (coords.empty()) {
        return {};
    }
    
    // Get MST
    auto mst = minimumSpanningTree(coords);
    if (mst.empty()) {
        return vector<int>(coords.size());
    }

    // Find vertices with odd degree
    map<int, int> degree;
    for (auto& e : mst) {
        if (e.first >= coords.size() || e.second >= coords.size()) {
            return vector<int>(coords.size());
        }
        degree[e.first]++;
        degree[e.second]++;
    }

    vector<int> odd;
    for (auto& [v, d] : degree) {
        if (d % 2 == 1) {
            if (v >= 0 && v < coords.size()) {
                odd.push_back(v);
            } else {
                return vector<int>(coords.size());
            }
        }
    }

    // Choose matching algorithm based on problem size
    vector<pair<int, int>> matching;
    if (odd.size() > 10000) {
        matching = greedyPerfectMatching(odd, coords);
    } else {
        // Create subgraph for Edmonds' algorithm
        vector<Point> odd_vertices;
        for (int idx : odd) {
            odd_vertices.push_back(coords[idx]);
        }
        
        auto edmonds_matching = blossom_edmonds(odd_vertices);
        
        // Convert Edmonds' matching indices back to original graph indices
        for (auto [u, v] : edmonds_matching) {
            if (u < odd.size() && v < odd.size()) {
                matching.emplace_back(odd[u], odd[v]);
            } else {
                return vector<int>(coords.size());
            }
        }
    }

    // Combine MST and matching edges into multigraph
    multimap<int, int> multigraph;
    try {
        for (auto& [u, v] : mst) {
            if (u >= coords.size() || v >= coords.size()) continue;
            multigraph.insert({u, v});
            multigraph.insert({v, u});
        }
        for (auto& [u, v] : matching) {
            if (u >= coords.size() || v >= coords.size()) continue;
            multigraph.insert({u, v});
            multigraph.insert({v, u});
        }
    } catch (const exception& e) {
        return vector<int>(coords.size());
    }

    if (multigraph.empty()) {
        return vector<int>(coords.size());
    }

    // Find Eulerian tour
    try {
        auto euler = eulerianTour(multigraph, 0);
        if (euler.empty()) {
            return vector<int>(coords.size());
        }
        
        // Make Hamiltonian
        auto hamilton = makeHamiltonian(euler);
        if (hamilton.empty()) {
            return vector<int>(coords.size());
        }
        
        // Verify path
        for (int v : hamilton) {
            if (v < 0 || v >= coords.size()) {
                return vector<int>(coords.size());
            }
        }
        
        return hamilton;
    } catch (const exception& e) {
        return vector<int>(coords.size());
    }
}

// Christofides Algorithm with Gabow's blossom matching
vector<int> christofidesPath_gabow(const vector<Point>& coords) {
    if (coords.empty()) {
        return {};
    }
    
    // Get MST
    auto mst = minimumSpanningTree(coords);
    if (mst.empty()) {
        return vector<int>(coords.size());
    }

    // Find vertices with odd degree
    map<int, int> degree;
    for (auto& e : mst) {
        if (e.first >= coords.size() || e.second >= coords.size()) {
            return vector<int>(coords.size());
        }
        degree[e.first]++;
        degree[e.second]++;
    }

    vector<int> odd;
    for (auto& [v, d] : degree) {
        if (d % 2 == 1) {
            if (v >= 0 && v < coords.size()) {
                odd.push_back(v);
            } else {
                return vector<int>(coords.size());
            }
        }
    }

    // Choose matching algorithm based on problem size
    vector<pair<int, int>> matching;
    if (odd.size() > 10000) {
        matching = greedyPerfectMatching(odd, coords);
    } else {
        // Create subgraph for Gabow's algorithm
        vector<Point> odd_vertices;
        for (int idx : odd) {
            odd_vertices.push_back(coords[idx]);
        }
        
        auto gabow_matching = blossom_gabow(odd_vertices);
        
        // Convert Gabow's matching indices back to original graph indices
        for (auto [u, v] : gabow_matching) {
            if (u < odd.size() && v < odd.size()) {
                matching.emplace_back(odd[u], odd[v]);
            } else {
                return vector<int>(coords.size());
            }
        }
    }

    // Combine MST and matching edges into multigraph
    multimap<int, int> multigraph;
    try {
        for (auto& [u, v] : mst) {
            if (u >= coords.size() || v >= coords.size()) continue;
            multigraph.insert({u, v});
            multigraph.insert({v, u});
        }
        for (auto& [u, v] : matching) {
            if (u >= coords.size() || v >= coords.size()) continue;
            multigraph.insert({u, v});
            multigraph.insert({v, u});
        }
    } catch (const exception& e) {
        return vector<int>(coords.size());
    }

    if (multigraph.empty()) {
        return vector<int>(coords.size());
    }

    // Find Eulerian tour
    try {
        auto euler = eulerianTour(multigraph, 0);
        if (euler.empty()) {
            return vector<int>(coords.size());
        }
        
        // Make Hamiltonian
        auto hamilton = makeHamiltonian(euler);
        if (hamilton.empty()) {
            return vector<int>(coords.size());
        }
        
        // Verify path
        for (int v : hamilton) {
            if (v < 0 || v >= coords.size()) {
                return vector<int>(coords.size());
            }
        }
        
        return hamilton;
    } catch (const exception& e) {
        return vector<int>(coords.size());
    }
}

