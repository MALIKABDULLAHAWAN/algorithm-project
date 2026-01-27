#include "graph.hpp"
#include <vector>
#include <queue>
#include <limits>

std::vector<CSRGraph::weight_t> dijkstra_seq(const CSRGraph& g, std::size_t source) {
    using T = CSRGraph::weight_t;
    const T INF = std::numeric_limits<T>::infinity();
    
    if (g.n == 0 || source >= g.n) return {};
    
    std::vector<T> dist(g.n, INF);
    dist[source] = 0;

    using Node = std::pair<T, std::size_t>; // (dist, node)
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;
    pq.emplace(0, source);

    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d != dist[u]) continue;
        for (std::size_t ei = g.row_offsets[u]; ei < g.row_offsets[u + 1]; ++ei) {
            std::size_t v = g.col_indices[ei];
            T w = g.weights[ei];
            if (w < 0) continue; // ignore negative weights for SSSP assumption
            if (dist[v] > d + w) {
                dist[v] = d + w;
                pq.emplace(dist[v], v);
            }
        }
    }
    return dist;
}
