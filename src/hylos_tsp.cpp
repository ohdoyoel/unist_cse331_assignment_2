#include "hylos_tsp.h"
#include <cmath>
#include <stack>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <random>

std::random_device rd;
std::mt19937 gen(rd());

HylosTSP::HylosTSP(const std::vector<Point>& cities, double alpha, bool debug)
    : cities_(cities), alpha_(alpha), debug_(debug) {}

std::vector<int> HylosTSP::solve() {
    debug_print("\n0. Starting Hylos TSP...\n");
    debug_print("1. Starting clustering...\n");
    cluster(ClusteringMethod::KMEANS);
    debug_print("   Created " + std::to_string(clusters_.size()) + " clusters\n");
    for (size_t i = 0; i < clusters_.size(); i++) {
        debug_print("   Cluster " + std::to_string(i) + " size: " + 
                   std::to_string(clusters_[i].cities.size()) + "\n");
    }
    
    debug_print("2. Determining cluster visit order...\n");
    cluster_order = determine_cluster_order();  // 멤버 변수에 저장
    debug_print("   Cluster visit order: ");
    for (int idx : cluster_order) {
        debug_print(std::to_string(idx) + " ");
    }
    debug_print("\n");
    
    debug_print("3. Solving cluster TSPs\n");
    solve_cluster_tsp();
    for (size_t i = 0; i < clusters_.size(); i++) {
        debug_print("   Cluster " + std::to_string(i) + " path size: " + 
                   std::to_string(clusters_[i].optimized_path.size()) + "\n");
    }
    
    debug_print("4. Finding connection points between clusters...\n");
    std::vector<int> final_path;
    // 각 클러스터의 연결점들을 찾고 경로 연결
    for (size_t i = 0; i < cluster_order.size(); i++) {
        int curr_idx = cluster_order[i];
        int prev_idx = cluster_order[(i + cluster_order.size() - 1) % cluster_order.size()];
        int next_idx = cluster_order[(i + 1) % cluster_order.size()];
        
        debug_print("    Finding connection points for cluster " + std::to_string(curr_idx) +
                   " (prev: " + std::to_string(prev_idx) +
                   ", next: " + std::to_string(next_idx) + ")\n");
        
        const auto& curr_cluster = clusters_[curr_idx];
        const auto& prev_cluster = clusters_[prev_idx];
        const auto& next_cluster = clusters_[next_idx];
        
        Point prev_centroid = get_cluster_centroid(prev_cluster);
        Point next_centroid = get_cluster_centroid(next_cluster);
        
        auto [entry_point, exit_point] = find_best_entry_exit_points(
            clusters_[curr_idx], prev_centroid, next_centroid);
            
        debug_print("   Selected entry point: " + std::to_string(entry_point) +
                   ", exit point: " + std::to_string(exit_point) + "\n");
        
        // 현재 클러스터의 reordered path를 final_path에 추가
        for (int city : clusters_[curr_idx].optimized_path) {
            final_path.push_back(city);
        }
    }

    debug_print("5. Final path: ");
    for (int city : final_path) {
        debug_print(std::to_string(city) + " ");
    }
    debug_print("\n");
    
    return final_path;
}

double HylosTSP::calculate_mean_std(const std::vector<double>& values, double& mean, double& std) {
    mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    double sq_sum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
    std = std::sqrt(sq_sum / values.size() - mean * mean);
    return std;
}

std::vector<int> HylosTSP::dfs_mst_order(const std::vector<std::pair<int, int>>& mst, const std::vector<Point>& points) {
    int n = points.size();  // 올바른 크기 사용
    std::vector<std::vector<int>> adj(n);
    std::vector<bool> visited(n, false);
    std::vector<int> order;
    
    // Build adjacency list from MST
    for (const auto& edge : mst) {
        adj[edge.first].push_back(edge.second);
        adj[edge.second].push_back(edge.first);
    }
    
    // DFS traversal
    std::function<void(int)> dfs = [&](int v) {
        visited[v] = true;
        order.push_back(v);
        for (int u : adj[v]) {
            if (!visited[u]) {
                dfs(u);
            }
        }
    };
    
    dfs(0);
    return order;
}

void HylosTSP::exclusive_clustering(const std::vector<int>& order) {
    std::vector<double> distances;
    for (size_t i = 1; i < order.size(); i++) {
        distances.push_back(euclideanDistance(cities_[order[i-1]], cities_[order[i]]));
    }
    
    double mean, std;
    calculate_mean_std(distances, mean, std);
    double threshold = mean + alpha_ * std;
    std::cerr << "   Distance threshold: " << threshold << "\n";
    
    // Initialize first cluster
    Cluster current_cluster;
    current_cluster.cities.push_back(order[0]);
    
    int min_cluster_size = 3;   // 최소 클러스터 크기
    int max_cluster_size = 22;  // 최대 클러스터 크기
    
    for (size_t i = 1; i < order.size(); i++) {
        double dist = euclideanDistance(cities_[order[i-1]], cities_[order[i]]);
        bool should_split = false;
        
        if (current_cluster.cities.size() >= max_cluster_size ||
            (dist > threshold && 
             current_cluster.cities.size() >= min_cluster_size && 
             order.size() - i >= min_cluster_size)) {
            should_split = true;
            std::cerr << "   Splitting at distance " << dist << " (threshold: " << threshold << ")\n";
        }
        
        if (should_split) {
            std::cerr << "   Created cluster with cities: ";
            for (int city : current_cluster.cities) {
                std::cerr << city << " ";
            }
            std::cerr << "\n";
            clusters_.push_back(current_cluster);
            current_cluster = Cluster();
        }
        
        current_cluster.cities.push_back(order[i]);
    }
    
    if (!current_cluster.cities.empty()) {
        if (current_cluster.cities.size() < min_cluster_size && !clusters_.empty()) {
            std::cerr << "   Merging last small cluster with previous cluster\n";
            auto& last_cluster = clusters_.back();
            last_cluster.cities.insert(last_cluster.cities.end(),
                                     current_cluster.cities.begin(),
                                     current_cluster.cities.end());
        } else {
            std::cerr << "   Created last cluster with cities: ";
            for (int city : current_cluster.cities) {
                std::cerr << city << " ";
            }
            std::cerr << "\n";
            clusters_.push_back(current_cluster);
        }
    }
}

void HylosTSP::solve_cluster_tsp() {
    for (size_t cluster_idx = 0; cluster_idx < clusters_.size(); cluster_idx++) {
        auto& cluster = clusters_[cluster_idx];
        
        std::cerr << "   Solving TSP for cluster " << cluster_idx 
                  << " (size: " << cluster.cities.size() << ")...\n";
        
        // 새로운 solve_cluster_tsp 함수를 사용하여 경로 계산
        cluster.optimized_path = solve_cluster_tsp(cluster);
        
        // 경로가 순환하도록 시작점을 마지막에 추가
        if (!cluster.optimized_path.empty()) {
            cluster.optimized_path.push_back(cluster.optimized_path[0]);
        }
    }
}

std::vector<std::pair<int, int>> HylosTSP::minimumSpanningTree(const std::vector<Point>& points) {
    int n = points.size();
    std::vector<std::pair<int, int>> mst;
    std::vector<bool> visited(n, false);
    std::vector<double> min_dist(n, std::numeric_limits<double>::max());
    std::vector<int> parent(n, -1);
    
    // Prim's algorithm
    min_dist[0] = 0;  // Start from vertex 0
    
    for (int i = 0; i < n; i++) {
        // Find vertex with minimum distance
        int u = -1;
        double curr_min = std::numeric_limits<double>::max();
        for (int j = 0; j < n; j++) {
            if (!visited[j] && min_dist[j] < curr_min) {
                curr_min = min_dist[j];
                u = j;
            }
        }
        
        if (u == -1) break;  // No more vertices to process
        visited[u] = true;
        
        // Add edge to MST
        if (parent[u] != -1) {
            mst.push_back({parent[u], u});
        }
        
        // Update distances
        for (int v = 0; v < n; v++) {
            if (!visited[v]) {
                double dist = euclideanDistance(points[u], points[v]);
                if (dist < min_dist[v]) {
                    min_dist[v] = dist;
                    parent[v] = u;
                }
            }
        }
    }
    
    return mst;
}

double HylosTSP::euclideanDistance(const Point& p1, const Point& p2) const {
    return std::sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

std::vector<Point> HylosTSP::initialize_centroids(int k) {
    std::vector<Point> centroids;
    // Use k-means++ initialization
    std::vector<double> distances(cities_.size(), std::numeric_limits<double>::max());
    
    // Choose first centroid randomly using modern C++ random
    std::uniform_int_distribution<int> dist(0, cities_.size() - 1);
    int first_idx = dist(gen);
    centroids.push_back(cities_[first_idx]);
    
    // Choose remaining centroids
    for (int i = 1; i < k; i++) {
        // Update distances
        for (size_t j = 0; j < cities_.size(); j++) {
            double min_dist = std::numeric_limits<double>::max();
            for (const auto& centroid : centroids) {
                min_dist = std::min(min_dist, euclideanDistance(cities_[j], centroid));
            }
            distances[j] = min_dist * min_dist;  // Square distances for probability
        }
        
        // Calculate sum for probability distribution
        double sum = std::accumulate(distances.begin(), distances.end(), 0.0);
        
        // Choose next centroid with probability proportional to distance using modern C++ random
        std::uniform_real_distribution<double> real_dist(0.0, sum);
        double r = real_dist(gen);
        double cumsum = 0.0;
        int chosen_idx = 0;
        for (size_t j = 0; j < cities_.size(); j++) {
            cumsum += distances[j];
            if (cumsum >= r) {
                chosen_idx = j;
                break;
            }
        }
        centroids.push_back(cities_[chosen_idx]);
    }
    return centroids;
}

std::vector<int> HylosTSP::assign_to_clusters(const std::vector<Point>& centroids) {
    std::vector<int> assignments(cities_.size());
    for (size_t i = 0; i < cities_.size(); i++) {
        double min_dist = std::numeric_limits<double>::max();
        int closest_centroid = 0;
        for (size_t j = 0; j < centroids.size(); j++) {
            double dist = euclideanDistance(cities_[i], centroids[j]);
            if (dist < min_dist) {
                min_dist = dist;
                closest_centroid = j;
            }
        }
        assignments[i] = closest_centroid;
    }
    return assignments;
}

std::vector<Point> HylosTSP::update_centroids(const std::vector<int>& assignments, int k) {
    std::vector<Point> new_centroids(k, {0.0, 0.0});
    std::vector<int> counts(k, 0);
    
    for (size_t i = 0; i < cities_.size(); i++) {
        int cluster = assignments[i];
        new_centroids[cluster].x += cities_[i].x;
        new_centroids[cluster].y += cities_[i].y;
        counts[cluster]++;
    }
    
    for (int i = 0; i < k; i++) {
        if (counts[i] > 0) {
            new_centroids[i].x /= counts[i];
            new_centroids[i].y /= counts[i];
        }
    }
    
    return new_centroids;
}

bool HylosTSP::has_converged(const std::vector<Point>& old_centroids, 
                            const std::vector<Point>& new_centroids,
                            double tolerance) {
    for (size_t i = 0; i < old_centroids.size(); i++) {
        if (euclideanDistance(old_centroids[i], new_centroids[i]) > tolerance) {
            return false;
        }
    }
    return true;
}

void HylosTSP::kmeans_clustering(const int k, const int MAX_SIZE, int max_iterations) {
    clusters_.clear();
    clusters_.resize(k);
    
    // Initialize centroids using k-means++
    std::vector<Point> centroids = initialize_centroids(k);
    std::vector<int> assignments(cities_.size(), -1);
    std::vector<std::vector<int>> cluster_cities(k);
    
    // Main k-means loop
    for (int iter = 0; iter < max_iterations; iter++) {
        bool changed = false;
        
        // Clear previous assignments
        for (auto& cluster : cluster_cities) {
            cluster.clear();
        }
        
        // Assign each city to nearest centroid
        for (size_t i = 0; i < cities_.size(); i++) {
            double min_dist = euclideanDistance(cities_[i], centroids[0]);
            int best_cluster = 0;
            
            for (int j = 1; j < k; j++) {
                double dist = euclideanDistance(cities_[i], centroids[j]);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_cluster = j;
                }
            }
            
            if (assignments[i] != best_cluster) {
                changed = true;
                assignments[i] = best_cluster;
            }
            cluster_cities[best_cluster].push_back(i);
        }
        
        if (!changed) break;
        
        // Update centroids
        for (int i = 0; i < k; i++) {
            if (cluster_cities[i].empty()) continue;
            
            Point new_centroid = {0, 0};
            for (int city : cluster_cities[i]) {
                new_centroid.x += cities_[city].x;
                new_centroid.y += cities_[city].y;
            }
            new_centroid.x /= cluster_cities[i].size();
            new_centroid.y /= cluster_cities[i].size();
            centroids[i] = new_centroid;
        }
    }
    
    // Split large clusters
    std::vector<Cluster> final_clusters;
    
    for (const auto& cities : cluster_cities) {
        if (cities.empty()) continue;
        
        if (cities.size() <= MAX_SIZE) {
            // Add cluster as is
            Cluster cluster;
            cluster.cities = cities;
            final_clusters.push_back(cluster);
        } else {
            // Split large cluster based on distances
            std::vector<std::pair<double, int>> distances;
            Point center = {0, 0};
            for (int city : cities) {
                center.x += cities_[city].x;
                center.y += cities_[city].y;
            }
            center.x /= cities.size();
            center.y /= cities.size();
            
            // Calculate distances from center
            for (int city : cities) {
                double dist = euclideanDistance(cities_[city], center);
                distances.push_back({dist, city});
            }
            
            // Sort by distance
            std::sort(distances.begin(), distances.end());
            
            // Create new clusters of size MAX_SIZE
            Cluster current_cluster;
            for (const auto& [dist, city] : distances) {
                current_cluster.cities.push_back(city);
                if (current_cluster.cities.size() == MAX_SIZE) {
                    final_clusters.push_back(current_cluster);
                    current_cluster = Cluster();
                }
            }
            
            // Add remaining cities if any
            if (!current_cluster.cities.empty()) {
                final_clusters.push_back(current_cluster);
            }
        }
    }
    
    clusters_ = final_clusters;
    
    // Verify no duplicates
    std::vector<bool> assigned(cities_.size(), false);
    for (const auto& cluster : clusters_) {
        for (int city : cluster.cities) {
            if (assigned[city]) {
                debug_print("Error: City " + std::to_string(city) + " is assigned to multiple clusters!\n");
            }
            assigned[city] = true;
        }
    }
    
    // Check for unassigned cities
    std::vector<int> unassigned;
    for (size_t i = 0; i < cities_.size(); i++) {
        if (!assigned[i]) {
            unassigned.push_back(i);
        }
    }
    
    // Add unassigned cities to nearest cluster that has room
    for (int city : unassigned) {
        double min_dist = std::numeric_limits<double>::max();
        size_t best_cluster = 0;
        
        for (size_t i = 0; i < clusters_.size(); i++) {
            if (clusters_[i].cities.size() >= MAX_SIZE) continue;
            
            Point center = {0, 0};
            for (int c : clusters_[i].cities) {
                center.x += cities_[c].x;
                center.y += cities_[c].y;
            }
            center.x /= clusters_[i].cities.size();
            center.y /= clusters_[i].cities.size();
            
            double dist = euclideanDistance(cities_[city], center);
            if (dist < min_dist) {
                min_dist = dist;
                best_cluster = i;
            }
        }
        
        clusters_[best_cluster].cities.push_back(city);
    }
}

Point HylosTSP::get_cluster_centroid(const Cluster& cluster) const {
    Point centroid = {0, 0};
    for (int city_idx : cluster.cities) {
        centroid.x += cities_[city_idx].x;
        centroid.y += cities_[city_idx].y;
    }
    centroid.x /= cluster.cities.size();
    centroid.y /= cluster.cities.size();
    return centroid;
}

std::vector<int> HylosTSP::determine_cluster_order() {
    std::vector<Point> centroids;
    for (const auto& cluster : clusters_) {
        centroids.push_back(get_cluster_centroid(cluster));
    }
    
    if (clusters_.size() <= 22) {
        // Use Held-Karp for small number of clusters
        std::cerr << "   Using Held-Karp for cluster ordering (size: " << clusters_.size() << ")\n";
        return heldkarpPath(centroids);
    } else {
        // Use Christofides for larger number of clusters
        std::cerr << "   Using Christofides for cluster ordering (size: " << clusters_.size() << ")\n";
        return christofidesPath(centroids);
    }
}

std::pair<int, int> HylosTSP::find_best_entry_exit_points(
    Cluster& cluster,
    const Point& prev_centroid,
    const Point& next_centroid) {
    
    int best_entry = cluster.optimized_path[0];
    int best_exit = cluster.optimized_path[1];
    double min_total_dist = std::numeric_limits<double>::max();

    // 마지막 점은 첫 점과 같으므로 제외하고 검사
    for (size_t i = 0; i < cluster.optimized_path.size() - 1; i++) {
        int curr = cluster.optimized_path[i];
        int next = cluster.optimized_path[(i + 1) % (cluster.optimized_path.size() - 1)];

        // 실제 도시 좌표를 사용하여 거리 계산
        double dist1 = euclideanDistance(cities_[curr], prev_centroid) +
                      euclideanDistance(cities_[next], next_centroid);
        
        double dist2 = euclideanDistance(cities_[next], prev_centroid) +
                      euclideanDistance(cities_[curr], next_centroid);

        if (dist1 < min_total_dist) {
            min_total_dist = dist1;
            best_entry = curr;
            best_exit = next;
        }
        if (dist2 < min_total_dist) {
            min_total_dist = dist2;
            best_entry = next;
            best_exit = curr;
        }
    }

    // 경로 재정렬
    std::vector<int> new_path;
    int entry_pos = -1, exit_pos = -1;
    
    // 실제 도시 인덱스로 위치 찾기
    for (size_t i = 0; i < cluster.optimized_path.size() - 1; i++) {
        if (cluster.optimized_path[i] == best_entry) entry_pos = i;
        if (cluster.optimized_path[i] == best_exit) exit_pos = i;
    }

    const int n = cluster.optimized_path.size() - 1;  // 마지막 점은 첫 점과 같으므로 제외

    bool entry_at_start = (entry_pos == 0);
    bool entry_at_end = (entry_pos == n-1);
    bool exit_at_start = (exit_pos == 0);
    bool exit_at_end = (exit_pos == n-1);

    // 경로 재구성 (기존 로직 유지)
    if (entry_at_start && exit_at_end) {
        for (int i = 0; i < n; i++) {
            new_path.push_back(cluster.optimized_path[i]);
        }
    } else if (exit_at_start && entry_at_end) {
        for (int i = n-1; i >= 0; i--) {
            new_path.push_back(cluster.optimized_path[i]);
        }
    } else if (entry_pos + 1 == exit_pos) {
        new_path.push_back(best_entry);
        
        for (int i = entry_pos - 1; i >= 0; i--) {
            new_path.push_back(cluster.optimized_path[i]);
        }
        
        for (int i = n-1; i > exit_pos; i--) {
            new_path.push_back(cluster.optimized_path[i]);
        }
        
        new_path.push_back(best_exit);
    } else if (exit_pos + 1 == entry_pos) {
        new_path.push_back(best_entry);
        
        for (int i = entry_pos + 1; i < n; i++) {
            new_path.push_back(cluster.optimized_path[i]);
        }
        
        for (int i = 0; i < exit_pos; i++) {
            new_path.push_back(cluster.optimized_path[i]);
        }
        
        new_path.push_back(best_exit);
    } else {
        new_path.push_back(best_entry);
        
        if (entry_pos < exit_pos) {
            for (int i = entry_pos + 1; i <= exit_pos; i++) {
                new_path.push_back(cluster.optimized_path[i]);
            }
        } else {
            for (int i = entry_pos - 1; i >= 0; i--) {
                new_path.push_back(cluster.optimized_path[i]);
            }
            for (int i = n - 1; i >= exit_pos; i--) {
                new_path.push_back(cluster.optimized_path[i]);
            }
        }
    }

    cluster.optimized_path = new_path;

    std::cerr << "   Reordered path: ";
    for (int city : cluster.optimized_path) {
        std::cerr << city << " ";
    }
    std::cerr << "\n";

    return {best_entry, best_exit};
}

std::vector<int> HylosTSP::mst2approx(const std::vector<Point>& points) {
    auto mst = minimumSpanningTree(points);
    return dfs_mst_order(mst, points);  // points 매개변수 전달
}

std::vector<int> HylosTSP::heldkarpPath(const std::vector<Point>& points) {
    int n = points.size();
    std::vector<std::vector<double>> dp(1 << n, std::vector<double>(n, 1e9));
    std::vector<std::vector<int>> next(1 << n, std::vector<int>(n, -1));
    
    dp[1][0] = 0;
    
    for (int mask = 0; mask < (1 << n); mask++) {
        for (int curr = 0; curr < n; curr++) {
            if (!(mask & (1 << curr))) continue;
            
            for (int next_city = 0; next_city < n; next_city++) {
                if (mask & (1 << next_city)) continue;
                
                int new_mask = mask | (1 << next_city);
                double new_cost = dp[mask][curr] + euclideanDistance(points[curr], points[next_city]);
                
                if (new_cost < dp[new_mask][next_city]) {
                    dp[new_mask][next_city] = new_cost;
                    next[new_mask][next_city] = curr;
                }
            }
        }
    }
    
    std::vector<int> path;
    int final_mask = (1 << n) - 1;
    int last_city = 0;
    double min_cost = std::numeric_limits<double>::max();
    
    for (int city = 0; city < n; city++) {
        double cost = dp[final_mask][city] + euclideanDistance(points[city], points[0]);
        if (cost < min_cost) {
            min_cost = cost;
            last_city = city;
        }
    }
    
    int curr_mask = final_mask;
    int curr_city = last_city;
    
    while (curr_city != -1) {
        path.push_back(curr_city);
        int prev_city = next[curr_mask][curr_city];
        if (prev_city == -1) break;
        curr_mask ^= (1 << curr_city);
        curr_city = prev_city;
    }
    
    std::reverse(path.begin(), path.end());
    return path;
}

DBSCANParams HylosTSP::find_optimal_eps_and_min_points() {
    const int n = cities_.size();
    std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n));
    
    // 1. 모든 도시 쌍 간의 거리 계산
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dist = euclideanDistance(cities_[i], cities_[j]);
            distance_matrix[i][j] = distance_matrix[j][i] = dist;
        }
    }
    
    // 2. 각 점에 대해 k-nearest neighbor 거리 계산 (k = 1 to 30)
    const int max_k = std::min(30, n-1);
    std::vector<std::vector<double>> k_distances(n, std::vector<double>(max_k));
    
    for (int i = 0; i < n; i++) {
        std::vector<double> distances = distance_matrix[i];
        std::sort(distances.begin(), distances.end());
        for (int k = 0; k < max_k; k++) {
            k_distances[i][k] = distances[k + 1];  // k+1 because first is always 0 (self)
        }
    }
    
    // 3. k-distance 그래프의 "elbow point" 찾기
    std::vector<double> avg_k_distances(max_k);
    for (int k = 0; k < max_k; k++) {
        double sum = 0;
        for (int i = 0; i < n; i++) {
            sum += k_distances[i][k];
        }
        avg_k_distances[k] = sum / n;
    }
    
    // 4. 기울기 변화가 가장 큰 지점 찾기 (elbow point)
    int optimal_k = 4;  // 기본값
    double max_angle_change = 0;
    
    for (int k = 1; k < max_k - 1; k++) {
        double prev_slope = avg_k_distances[k] - avg_k_distances[k-1];
        double next_slope = avg_k_distances[k+1] - avg_k_distances[k];
        double angle_change = std::abs(std::atan(next_slope) - std::atan(prev_slope));
        
        if (angle_change > max_angle_change) {
            max_angle_change = angle_change;
            optimal_k = k;
        }
    }
    
    // 5. optimal_k를 사용하여 eps 계산
    double eps = 0;
    for (int i = 0; i < n; i++) {
        eps += k_distances[i][optimal_k];
    }
    eps /= n;
    
    // 6. min_points 결정
    // min_points는 optimal_k + 1로 설정 (자기 자신 포함)
    int min_points = optimal_k + 1;
    
    // 7. 제한 조건 적용
    min_points = std::min(std::max(min_points, 4), 20);  // 4 ≤ min_points ≤ 20
    
    std::cerr << "   Optimal parameters found:\n"
              << "     k-distance elbow point: " << optimal_k << "\n"
              << "     eps: " << eps << "\n"
              << "     min_points: " << min_points << "\n";
    
    return {eps, min_points};
}

void HylosTSP::cluster(ClusteringMethod method, double eps, int min_points) {
    switch (method) {
        case ClusteringMethod::KMEANS: {
            // n이 작을 때는 더 작은 클러스터를, n이 클 때는 더 큰 클러스터를 만듭니다.
            // sqrt(n)을 기준으로 MAX_SIZE를 조정합니다.
            const int n = cities_.size();
            const int sqrt_n = static_cast<int>(std::sqrt(n));
            
            // MAX_SIZE는 sqrt(n)의 2~3배 정도로 설정
            // 단, 최소 22으로 제한
            const int MAX_SIZE = std::max(22, 2 * sqrt_n);
            
            // k는 n/MAX_SIZE를 올림한 값으로 설정
            const int k = (n + MAX_SIZE - 1) / MAX_SIZE;  // ceiling division
            
            std::cerr << "   Total cities: " << n << "\n"
                      << "   MAX_SIZE: " << MAX_SIZE << "\n"
                      << "   Number of clusters (k): " << k << "\n";
            
            kmeans_clustering(k, MAX_SIZE);
            break;
        }
        case ClusteringMethod::DBSCAN: {
            auto params = find_optimal_eps_and_min_points();
            std::cerr << "   Using DBSCAN with:\n"
                      << "     eps=" << params.eps << "\n"
                      << "     min_points=" << params.min_points << "\n";
            dbscan_clustering(params.eps, params.min_points);
            break;
        }
    }
}

void HylosTSP::dbscan_clustering(double eps, int min_points) {
    std::cerr << "Starting DBSCAN clustering with eps=" << eps 
              << ", min_points=" << min_points << "...\n";
              
    clusters_.clear();
    std::vector<bool> visited(cities_.size(), false);
    std::vector<int> point_labels(cities_.size(), -1);  // -1: unassigned, -2: noise, >=0: cluster index
    
    int cluster_idx = 0;
    for (size_t i = 0; i < cities_.size(); i++) {
        if (visited[i] || point_labels[i] != -1) continue;
        
        visited[i] = true;
        std::vector<int> neighbors = region_query(i, eps);
        
        if (neighbors.size() < min_points) {
            // Mark as noise
            point_labels[i] = -2;
            continue;
        }
        
        // Start a new cluster
        Cluster new_cluster;
        new_cluster.cities.push_back(i);
        point_labels[i] = cluster_idx;
        
        // Process neighbors
        std::queue<int> queue;
        for (int neighbor : neighbors) {
            if (!visited[neighbor]) {
                queue.push(neighbor);
            }
        }
        
        while (!queue.empty()) {
            int current = queue.front();
            queue.pop();
            
            if (!visited[current]) {
                visited[current] = true;
                std::vector<int> current_neighbors = region_query(current, eps);
                
                if (current_neighbors.size() >= min_points) {
                    for (int neighbor : current_neighbors) {
                        if (!visited[neighbor]) {
                            queue.push(neighbor);
                        }
                    }
                }
            }
            
            // Add to cluster if not already in another cluster
            if (point_labels[current] == -1) {
                new_cluster.cities.push_back(current);
                point_labels[current] = cluster_idx;
            }
        }
        
        clusters_.push_back(new_cluster);
        cluster_idx++;
    }
    
    // Handle noise points: create singleton clusters
    for (size_t i = 0; i < cities_.size(); i++) {
        if (point_labels[i] == -2 || point_labels[i] == -1) {
            Cluster noise_cluster;
            noise_cluster.cities.push_back(i);
            clusters_.push_back(noise_cluster);
        }
    }
    
    // Print cluster information
    std::cerr << "Created " << clusters_.size() << " clusters:\n";
    for (size_t i = 0; i < clusters_.size(); i++) {
        std::cerr << "Cluster " << i << " size: " << clusters_[i].cities.size() << "\n";
    }
}

std::vector<int> HylosTSP::region_query(int point_idx, double eps) const {
    std::vector<int> neighbors;
    for (size_t i = 0; i < cities_.size(); i++) {
        if (i != point_idx && 
            euclideanDistance(cities_[point_idx], cities_[i]) <= eps) {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

std::vector<int> HylosTSP::solve_cluster_tsp(const Cluster& cluster) {
    std::vector<Point> cluster_points;
    // 원래 도시 인덱스를 저장
    std::vector<int> original_indices = cluster.cities;
    
    for (int city_idx : cluster.cities) {
        cluster_points.push_back(cities_[city_idx]);
    }
    
    std::vector<int> local_path;
    if (cluster.cities.size() <= 22) {
        std::cerr << "   Using Held-Karp for cluster size " << cluster.cities.size() << "\n";
        local_path = heldkarpPath(cluster_points);
    } else {
        std::cerr << "   Using Christofides for cluster size " << cluster.cities.size() << "\n";
        local_path = christofidesPath(cluster_points);
    }
    
    // 로컬 인덱스를 원래 도시 인덱스로 변환
    std::vector<int> global_path;
    for (int local_idx : local_path) {
        global_path.push_back(original_indices[local_idx]);
    }
    
    std::cerr << "   Path: ";
    for (int city : global_path) {
        std::cerr << city << " ";
    }
    std::cerr << "\n";
    
    return global_path;
} 