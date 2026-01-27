#include "graph.hpp"
#include <vector>
#include <limits>
#include <cmath>
#include <queue>
#include <utility>
#include <algorithm>
#ifdef HAS_OPENMP
#include <omp.h>
#endif

// Parallel Δ-stepping (Meyer & Sanders) for non-negative weights.
// Buckets distances into ranges [k*delta, (k+1)*delta). Processes light edges inside a bucket
// until closure, then relaxes heavy edges to future buckets.
std::vector<CSRGraph::weight_t> delta_stepping_omp(const CSRGraph& g, std::size_t source, double delta, int threads) {
    using T = CSRGraph::weight_t;
    const T INF = std::numeric_limits<T>::infinity();

    if (g.n == 0 || source >= g.n) return {};
    if (delta <= 0) delta = 1.0; // fallback

    std::vector<T> dist(g.n, INF);
    dist[source] = 0.0;

    auto bucket_index = [&](T d) -> std::size_t {
        if (d == INF) return static_cast<std::size_t>(-1);
        return static_cast<std::size_t>(std::floor(d / delta));
    };

    // Dynamic bucket structure: vector of vectors, grown as needed
    std::vector<std::vector<std::size_t>> B(1);
    B[0].push_back(source);

    std::size_t i = 0; // current bucket index

#ifdef HAS_OPENMP
    if (threads > 0) omp_set_num_threads(threads);
#endif

    auto relax_edge = [&](std::size_t u, std::size_t ei, std::vector<std::pair<std::size_t,T>>& relaxes) {
        std::size_t v = g.col_indices[ei];
        T w = g.weights[ei];
        if (w < 0) return; // skip negative weights (algorithm assumes non-negative)
        T nd = dist[u] + w;
        if (nd < dist[v]) {
            relaxes.emplace_back(v, nd);
        }
    };

    // Process buckets in increasing order
    while (i < B.size()) {
        // Skip empty buckets
        if (B[i].empty()) { ++i; continue; }

        // R is the set of vertices whose light edges we still need to process at this bucket
        std::vector<std::size_t> R;
        R.swap(B[i]);

        // S accumulates vertices settled for this bucket (after light-edge closure)
        std::vector<std::size_t> S;

        // Closure over light edges
        while (!R.empty()) {
            // Move R into S
            S.insert(S.end(), R.begin(), R.end());

            // Gather relaxations from light edges in parallel
            std::vector<std::pair<std::size_t,T>> relaxes; // (v, new_dist)

#ifdef HAS_OPENMP
            #pragma omp parallel
            {
                std::vector<std::pair<std::size_t,T>> local_relaxes;
                #pragma omp for schedule(dynamic, 1024) nowait
                for (long long idx = 0; idx < static_cast<long long>(R.size()); ++idx) {
                    std::size_t u = R[static_cast<std::size_t>(idx)];
                    for (std::size_t ei = g.row_offsets[u]; ei < g.row_offsets[u + 1]; ++ei) {
                        if (g.weights[ei] <= delta) {
                            relax_edge(u, ei, local_relaxes);
                        }
                    }
                }
                #pragma omp critical
                relaxes.insert(relaxes.end(), local_relaxes.begin(), local_relaxes.end());
            }
#else
            for (std::size_t u : R) {
                for (std::size_t ei = g.row_offsets[u]; ei < g.row_offsets[u + 1]; ++ei) {
                    if (g.weights[ei] <= delta) {
                        relax_edge(u, ei, relaxes);
                    }
                }
            }
#endif

            // Apply relaxations with atomic-like effect by comparing and updating distances
            // Also build the next R' consisting of vertices still in current bucket after update
            std::vector<std::size_t> nextR;

            // Sort by vertex to reduce contention and keep best nd
            std::sort(relaxes.begin(), relaxes.end());
            std::size_t k = 0;
            while (k < relaxes.size()) {
                std::size_t v = relaxes[k].first;
                T best = relaxes[k].second;
                std::size_t j = k + 1;
                while (j < relaxes.size() && relaxes[j].first == v) {
                    best = std::min(best, relaxes[j].second);
                    ++j;
                }
                if (best < dist[v]) {
                    dist[v] = best;
                    std::size_t bi = bucket_index(best);
                    if (bi == i) {
                        nextR.push_back(v);
                    } else {
                        // Light edge relaxation to a future bucket
                        if (bi >= B.size()) B.resize(bi + 1);
                        B[bi].push_back(v);
                    }
                }
                k = j;
            }
            R.swap(nextR);
        }

        // Heavy edges: relax from S to future buckets
        std::vector<std::pair<std::size_t,T>> heavy_relaxes;

#ifdef HAS_OPENMP
        #pragma omp parallel
        {
            std::vector<std::pair<std::size_t,T>> local_relaxes;
            #pragma omp for schedule(dynamic, 1024) nowait
            for (long long idx = 0; idx < static_cast<long long>(S.size()); ++idx) {
                std::size_t u = S[static_cast<std::size_t>(idx)];
                for (std::size_t ei = g.row_offsets[u]; ei < g.row_offsets[u + 1]; ++ei) {
                    if (g.weights[ei] > delta) {
                        relax_edge(u, ei, local_relaxes);
                    }
                }
            }
            #pragma omp critical
            heavy_relaxes.insert(heavy_relaxes.end(), local_relaxes.begin(), local_relaxes.end());
        }
#else
        for (std::size_t u : S) {
            for (std::size_t ei = g.row_offsets[u]; ei < g.row_offsets[u + 1]; ++ei) {
                if (g.weights[ei] > delta) {
                    relax_edge(u, ei, heavy_relaxes);
                }
            }
        }
#endif

        // Coalesce heavy relaxations per vertex (keep best)
        std::sort(heavy_relaxes.begin(), heavy_relaxes.end());
        std::size_t k2 = 0;
        while (k2 < heavy_relaxes.size()) {
            std::size_t v = heavy_relaxes[k2].first;
            T best = heavy_relaxes[k2].second;
            std::size_t j2 = k2 + 1;
            while (j2 < heavy_relaxes.size() && heavy_relaxes[j2].first == v) {
                best = std::min(best, heavy_relaxes[j2].second);
                ++j2;
            }
            if (best < dist[v]) {
                // Update dist and enqueue into appropriate future bucket
                dist[v] = best;
                std::size_t bi = bucket_index(best);
                if (bi >= B.size()) B.resize(bi + 1);
                B[bi].push_back(v);
            }
            k2 = j2;
        }

        // Move to next bucket
        ++i;
    }

    return dist;
}
