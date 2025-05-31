#include "hylos_tsp.h"
#include <cmath>
#include <stack>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <functional>
#include <iostream>
#include <limits>

HylosTSP::HylosTSP(const std::vector<Point>& cities, double alpha)
    : cities_(cities), alpha_(alpha) {}

std::vector<int> HylosTSP::solve() {
    std::cout << "1. Starting k-means clustering...\n";
    int k = cities_.size() / 15;  // Aim for ~15 cities per cluster
    k = std::max(10, std::min(k, 20));  // Keep k between 10 and 20
    kmeans_clustering(k);
    std::cout << "   Created " << clusters_.size() << " clusters\n";
    for (size_t i = 0; i < clusters_.size(); i++) {
        std::cout << "   Cluster " << i << " size: " << clusters_[i].cities.size() << "\n";
    }
    
    std::cout << "2. Determining cluster visit order...\n";
    cluster_order = determine_cluster_order();  // 멤버 변수에 저장
    std::cout << "   Cluster visit order: ";
    for (int idx : cluster_order) {
        std::cout << idx << " ";
    }
    std::cout << "\n";
    
    std::cout << "3. Solving cluster TSPs using Held-Karp...\n";
    solve_cluster_tsp();
    for (size_t i = 0; i < clusters_.size(); i++) {
        std::cout << "   Cluster " << i << " path size: " << clusters_[i].optimized_path.size() << "\n";
    }
    
    std::cout << "4. Finding connection points between clusters...\n";
    std::vector<int> final_path;
    
    // 각 클러스터의 연결점들을 찾고 경로 연결
    for (size_t i = 0; i < cluster_order.size(); i++) {
        int curr_idx = cluster_order[i];
        int prev_idx = cluster_order[(i + cluster_order.size() - 1) % cluster_order.size()];
        int next_idx = cluster_order[(i + 1) % cluster_order.size()];
        
        std::cout << "\nFinding connection points for cluster " << curr_idx 
                  << " (prev: " << prev_idx 
                  << ", next: " << next_idx << ")\n";
        
        const auto& curr_cluster = clusters_[curr_idx];
        const auto& prev_cluster = clusters_[prev_idx];
        const auto& next_cluster = clusters_[next_idx];
        
        Point prev_centroid = get_cluster_centroid(prev_cluster);
        Point next_centroid = get_cluster_centroid(next_cluster);
        
        auto [entry_point, exit_point] = find_best_entry_exit_points(
            clusters_[curr_idx], prev_centroid, next_centroid);
            
        std::cout << "   Selected entry point: " << entry_point 
                  << ", exit point: " << exit_point << "\n";
        
        // 현재 클러스터의 reordered path를 final_path에 추가
        for (int city : clusters_[curr_idx].optimized_path) {
            final_path.push_back(city);
        }
    }

    std::cout << "\n5. Final path: ";
    for (int city : final_path) {
        std::cout << city << " ";
    }
    std::cout << "\n";
    
    // // 순환 경로 완성을 위해 첫 번째 점 추가
    // final_path.push_back(final_path[0]);
    
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
    std::cout << "   Distance threshold: " << threshold << "\n";
    
    // Initialize first cluster
    Cluster current_cluster;
    current_cluster.cities.push_back(order[0]);
    
    int min_cluster_size = 3;   // 최소 클러스터 크기
    int max_cluster_size = 20;  // 최대 클러스터 크기
    
    for (size_t i = 1; i < order.size(); i++) {
        double dist = euclideanDistance(cities_[order[i-1]], cities_[order[i]]);
        bool should_split = false;
        
        if (current_cluster.cities.size() >= max_cluster_size ||
            (dist > threshold && 
             current_cluster.cities.size() >= min_cluster_size && 
             order.size() - i >= min_cluster_size)) {
            should_split = true;
            std::cout << "   Splitting at distance " << dist << " (threshold: " << threshold << ")\n";
        }
        
        if (should_split) {
            std::cout << "   Created cluster with cities: ";
            for (int city : current_cluster.cities) {
                std::cout << city << " ";
            }
            std::cout << "\n";
            clusters_.push_back(current_cluster);
            current_cluster = Cluster();
        }
        
        current_cluster.cities.push_back(order[i]);
    }
    
    if (!current_cluster.cities.empty()) {
        if (current_cluster.cities.size() < min_cluster_size && !clusters_.empty()) {
            std::cout << "   Merging last small cluster with previous cluster\n";
            auto& last_cluster = clusters_.back();
            last_cluster.cities.insert(last_cluster.cities.end(),
                                     current_cluster.cities.begin(),
                                     current_cluster.cities.end());
        } else {
            std::cout << "   Created last cluster with cities: ";
            for (int city : current_cluster.cities) {
                std::cout << city << " ";
            }
            std::cout << "\n";
            clusters_.push_back(current_cluster);
        }
    }
}

void HylosTSP::solve_cluster_tsp() {
    for (size_t cluster_idx = 0; cluster_idx < clusters_.size(); cluster_idx++) {
        auto& cluster = clusters_[cluster_idx];
        int n = cluster.cities.size();
        std::cout << "   Solving TSP for cluster " << cluster_idx << " of size " << n << "...\n";
        
        if (n == 1) {
            // Single city case
            cluster.optimized_path = {cluster.cities[0]};
        } else if (n == 2) {
            // Two cities case
            cluster.optimized_path = {cluster.cities[0], cluster.cities[1], cluster.cities[0]};
        } else {
            // Held-Karp DP를 사용한 최적 경로 찾기
            std::vector<std::vector<double>> dp(1 << n, std::vector<double>(n, 1e9));
            std::vector<std::vector<int>> next(1 << n, std::vector<int>(n, -1));
            
            dp[1][0] = 0;
            
            for (int mask = 0; mask < (1 << n); mask++) {
                for (int curr = 0; curr < n; curr++) {
                    if (!(mask & (1 << curr))) continue;
                    
                    for (int next_city = 0; next_city < n; next_city++) {
                        if (mask & (1 << next_city)) continue;
                        
                        int new_mask = mask | (1 << next_city);
                        double new_cost = dp[mask][curr] + 
                                        euclideanDistance(cities_[cluster.cities[curr]], 
                                                        cities_[cluster.cities[next_city]]);
                        
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
                double cost = dp[final_mask][city] + 
                             euclideanDistance(cities_[cluster.cities[city]], 
                                             cities_[cluster.cities[0]]);
                if (cost < min_cost) {
                    min_cost = cost;
                    last_city = city;
                }
            }
            
            int curr_mask = final_mask;
            int curr_city = last_city;
            
            while (curr_city != -1) {
                path.push_back(cluster.cities[curr_city]);
                int prev_city = next[curr_mask][curr_city];
                if (prev_city == -1) break;
                curr_mask ^= (1 << curr_city);
                curr_city = prev_city;
            }
            
            std::reverse(path.begin(), path.end());
            path.push_back(path[0]);  // 순환 경로 완성
            
            cluster.optimized_path = path;
        }
        
        std::cout << "   Cluster " << cluster_idx << " optimized path: ";
        for (int city : cluster.optimized_path) {
            std::cout << city << " ";
        }
        std::cout << "\n";
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
    
    // Choose first centroid randomly
    int first_idx = rand() % cities_.size();
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
        
        // Choose next centroid with probability proportional to distance
        double r = (double)rand() / RAND_MAX * sum;
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

void HylosTSP::kmeans_clustering(int k, int max_iterations) {
    clusters_.clear();
    clusters_.resize(k);
    
    // Initialize centroids by selecting k evenly spaced cities
    std::vector<Point> centroids;
    int step = cities_.size() / k;
    for (int i = 0; i < k; i++) {
        centroids.push_back(cities_[i * step]);
    }
    
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
    
    // Split large clusters and merge small ones
    const int MAX_SIZE = 20;
    const int MIN_SIZE = 2;
    std::vector<std::vector<int>> new_clusters;
    
    for (const auto& cluster : cluster_cities) {
        if (cluster.size() < MIN_SIZE) continue;  // Will handle small clusters later
        
        if (cluster.size() <= MAX_SIZE) {
            new_clusters.push_back(cluster);
        } else {
            // Split large cluster based on distance from centroid
            std::vector<std::pair<double, int>> distances;
            Point center = {0, 0};
            for (int city : cluster) {
                center.x += cities_[city].x;
                center.y += cities_[city].y;
            }
            center.x /= cluster.size();
            center.y /= cluster.size();
            
            for (int city : cluster) {
                double dist = euclideanDistance(cities_[city], center);
                distances.push_back({dist, city});
            }
            
            // Sort by distance
            std::sort(distances.begin(), distances.end());
            
            // Create new clusters of size MAX_SIZE
            std::vector<int> current_cluster;
            for (const auto& [dist, city] : distances) {
                current_cluster.push_back(city);
                if (current_cluster.size() == MAX_SIZE) {
                    new_clusters.push_back(current_cluster);
                    current_cluster.clear();
                }
            }
            
            // Add remaining cities if any
            if (!current_cluster.empty()) {
                new_clusters.push_back(current_cluster);
            }
        }
    }
    
    // Handle small clusters by merging them with the nearest cluster
    for (const auto& cluster : cluster_cities) {
        if (cluster.size() < MIN_SIZE) {
            for (int city : cluster) {
                double min_dist = std::numeric_limits<double>::max();
                size_t best_cluster = 0;
                
                // Find nearest cluster that has room
                for (size_t i = 0; i < new_clusters.size(); i++) {
                    if (new_clusters[i].size() >= MAX_SIZE) continue;
                    
                    Point center = {0, 0};
                    for (int c : new_clusters[i]) {
                        center.x += cities_[c].x;
                        center.y += cities_[c].y;
                    }
                    center.x /= new_clusters[i].size();
                    center.y /= new_clusters[i].size();
                    
                    double dist = euclideanDistance(cities_[city], center);
                    if (dist < min_dist) {
                        min_dist = dist;
                        best_cluster = i;
                    }
                }
                
                if (new_clusters[best_cluster].size() < MAX_SIZE) {
                    new_clusters[best_cluster].push_back(city);
                }
            }
        }
    }
    
    // Create final clusters
    clusters_.clear();
    for (const auto& cities : new_clusters) {
        if (!cities.empty()) {
            Cluster cluster;
            cluster.cities = cities;
            clusters_.push_back(cluster);
        }
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
    
    if (clusters_.size() <= 20) {
        // Use Held-Karp for small number of clusters
        std::cout << "   Using Held-Karp for cluster ordering\n";
        return heldkarpPath(centroids);
    } else {
        // Use MST 2-approximation for larger number of clusters
        std::cout << "   Using MST 2-approximation for cluster ordering\n";
        return mst2approx(centroids);
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

    int entry_pos = -1, exit_pos = -1;
    for (size_t i = 0; i < cluster.optimized_path.size() - 1; i++) {
        if (cluster.optimized_path[i] == best_entry) entry_pos = i;
        if (cluster.optimized_path[i] == best_exit) exit_pos = i;
    }

    std::vector<int> new_path;
    const int n = cluster.optimized_path.size() - 1;  // 마지막 점은 첫 점과 같으므로 제외

    bool entry_at_start = (entry_pos == 0);
    bool entry_at_end = (entry_pos == n-1);
    bool exit_at_start = (exit_pos == 0);
    bool exit_at_end = (exit_pos == n-1);

    // Case 1: a, ..., b (entry가 시작점이고 exit가 끝점인 경우)
    if (entry_at_start && exit_at_end) {
        for (int i = 0; i < n; i++) {
            new_path.push_back(cluster.optimized_path[i]);
        }
    }
    // Case 2: b, ..., a (exit가 시작점이고 entry가 끝점인 경우)
    else if (exit_at_start && entry_at_end) {
        for (int i = n-1; i >= 0; i--) {
            new_path.push_back(cluster.optimized_path[i]);
        }
    }
    // Case 3: ...1, a, b, ...2 (entry와 exit가 연속)
    else if (entry_pos + 1 == exit_pos) {
        new_path.push_back(best_entry);  // a
        
        // reverse(...1)
        for (int i = entry_pos - 1; i >= 0; i--) {
            new_path.push_back(cluster.optimized_path[i]);
        }
        
        // reverse(...2)
        for (int i = n-1; i > exit_pos; i--) {
            new_path.push_back(cluster.optimized_path[i]);
        }
        
        new_path.push_back(best_exit);  // b
    }
    // Case 4: ...1, b, a, ...2 (exit와 entry가 연속)
    else if (exit_pos + 1 == entry_pos) {
        new_path.push_back(best_entry);  // a
        
        // ...2
        for (int i = entry_pos + 1; i < n; i++) {
            new_path.push_back(cluster.optimized_path[i]);
        }
        
        // ...1
        for (int i = 0; i < exit_pos; i++) {
            new_path.push_back(cluster.optimized_path[i]);
        }
        
        new_path.push_back(best_exit);  // b
    }
    // 일반적인 경우: 가까운 쪽으로 이동
    else {
        new_path.push_back(best_entry);
        
        if (entry_pos < exit_pos) {
            // entry -> ... -> exit 순서로 있는 경우
            for (int i = entry_pos + 1; i <= exit_pos; i++) {
                new_path.push_back(cluster.optimized_path[i]);
            }
        } else {
            // exit -> ... -> entry 순서로 있는 경우
            for (int i = entry_pos - 1; i >= 0; i--) {
                new_path.push_back(cluster.optimized_path[i]);
            }
            for (int i = n - 1; i >= exit_pos; i--) {
                new_path.push_back(cluster.optimized_path[i]);
            }
        }
    }

    cluster.optimized_path = new_path;

    std::cout << "   Reordered path: ";
    for (int city : cluster.optimized_path) {
        std::cout << city << " ";
    }
    std::cout << "\n";

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