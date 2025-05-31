#pragma once

#include "tsp.h"
#include <vector>
#include <utility>
#include <cmath>

struct Cluster {
    std::vector<int> cities;
    std::vector<int> optimized_path;
};

class HylosTSP {
public:
    HylosTSP(const std::vector<Point>& cities, double alpha = 1.0);
    std::vector<int> solve();
    const std::vector<Cluster>& get_clusters() const { return clusters_; }
    Point get_cluster_centroid(const Cluster& cluster) const;

private:
    std::vector<Point> cities_;
    std::vector<Cluster> clusters_;
    std::vector<int> cluster_order;  // 클러스터 방문 순서
    double alpha_;

    // MST와 DFS 관련 함수들
    std::vector<std::pair<int, int>> minimumSpanningTree(const std::vector<Point>& points);
    std::vector<int> dfs_mst_order(const std::vector<std::pair<int, int>>& mst, const std::vector<Point>& points);
    std::vector<int> mst2approx(const std::vector<Point>& points);
    
    // 클러스터링 및 TSP 해결 함수들
    void exclusive_clustering(const std::vector<int>& order);
    void solve_cluster_tsp();
    std::vector<int> heldkarpPath(const std::vector<Point>& points);
    
    // 유틸리티 함수들
    double calculate_mean_std(const std::vector<double>& values, double& mean, double& std);
    double euclideanDistance(const Point& p1, const Point& p2) const;
    
    // k-means clustering methods
    void kmeans_clustering(int k, int max_iterations = 100);
    std::vector<Point> initialize_centroids(int k);
    std::vector<int> assign_to_clusters(const std::vector<Point>& centroids);
    std::vector<Point> update_centroids(const std::vector<int>& assignments, int k);
    bool has_converged(const std::vector<Point>& old_centroids, 
                      const std::vector<Point>& new_centroids,
                      double tolerance = 1e-4);
                      
    // Cluster handling methods
    std::vector<int> determine_cluster_order();
    std::pair<int, int> find_best_entry_exit_points(Cluster& cluster, const Point& prev_centroid, const Point& next_centroid);
}; 