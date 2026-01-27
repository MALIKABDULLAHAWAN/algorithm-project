#include "graph.hpp"
#include <algorithm>
#include <tuple>
#include <unordered_set>

static inline double urand(std::mt19937& rng) {
    static thread_local std::uniform_real_distribution<double> U(0.0, 1.0);
    return U(rng);
}

CSRGraph build_csr(std::size_t n, const std::vector<std::tuple<std::size_t,std::size_t,double>>& edges) {
    CSRGraph g;
    g.n = n;
    g.m = edges.size();
    g.row_offsets.assign(n + 1, 0);
    g.col_indices.resize(g.m);
    g.weights.resize(g.m);

    for (const auto& e : edges) {
        std::size_t u = std::get<0>(e);
        if (u >= n) continue; // ignore invalid
        g.row_offsets[u + 1]++;
    }
    for (std::size_t i = 1; i <= n; ++i) g.row_offsets[i] += g.row_offsets[i - 1];

    std::vector<std::size_t> cur = g.row_offsets; // cursor per row
    for (const auto& e : edges) {
        std::size_t u = std::get<0>(e);
        std::size_t v = std::get<1>(e);
        double w = std::get<2>(e);
        if (u >= n || v >= n) continue; // ignore invalid
        auto idx = cur[u]++;
        g.col_indices[idx] = v;
        g.weights[idx] = w;
    }
    return g;
}

CSRGraph generate_er(std::size_t n, double p, unsigned int seed, double wmin, double wmax, bool undirected) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> pr(0.0, 1.0);
    std::uniform_real_distribution<double> wr(wmin, wmax);

    std::vector<std::tuple<std::size_t,std::size_t,double>> edges;
    edges.reserve(static_cast<std::size_t>(p * n * (n - 1)) + 1);

    for (std::size_t u = 0; u < n; ++u) {
        for (std::size_t v = 0; v < n; ++v) {
            if (u == v) continue;
            if (pr(rng) < p) {
                double w = wr(rng);
                edges.emplace_back(u, v, w);
            }
        }
    }

    if (undirected) {
        // ensure symmetry by duplicating reverse edges when missing (approximate)
        // For simplicity, we just add reverse edges for each added edge.
        std::vector<std::tuple<std::size_t,std::size_t,double>> und;
        und.reserve(edges.size() * 2);
        for (auto &e : edges) {
            auto u = std::get<0>(e); auto v = std::get<1>(e); auto w = std::get<2>(e);
            und.emplace_back(u, v, w);
            und.emplace_back(v, u, w);
        }
        edges.swap(und);
    }

    return build_csr(n, edges);
}

CSRGraph generate_ba(std::size_t n, std::size_t m0, std::size_t m_attach, unsigned int seed,
                     double wmin, double wmax, bool undirected) {
    if (n == 0) return {};
    if (m0 < 1) m0 = 1;
    if (m_attach < 1) m_attach = 1;
    if (m0 >= n) m0 = std::min<std::size_t>(n, std::max<std::size_t>(1, m0));

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> wr(wmin, wmax);

    std::vector<std::tuple<std::size_t,std::size_t,double>> edges;
    edges.reserve(n * m_attach * (undirected ? 2 : 1));

    // Start with m0 fully connected (clique)
    for (std::size_t u = 0; u < m0; ++u) {
        for (std::size_t v = u + 1; v < m0; ++v) {
            double w = wr(rng);
            edges.emplace_back(u, v, w);
            edges.emplace_back(v, u, w); // always add both directions for clique
        }
    }

    // Preferential attachment: maintain list of endpoints proportional to degree
    std::vector<std::size_t> targets;
    targets.reserve(m0 * m0 * 2 + n * m_attach * 2);
    // initialize targets with clique degrees
    for (std::size_t u = 0; u < m0; ++u) {
        for (std::size_t k = 0; k < (m0 - 1) * (undirected ? 2 : 1); ++k) targets.push_back(u);
    }

    for (std::size_t u = m0; u < n; ++u) {
        // choose m_attach unique neighbors based on targets
        std::unordered_set<std::size_t> chosen;
        while (chosen.size() < m_attach && !targets.empty()) {
            std::size_t idx = static_cast<std::size_t>(urand(rng) * targets.size());
            if (idx >= targets.size()) idx = targets.size() - 1;
            std::size_t v = targets[idx];
            if (v == u) continue;
            chosen.insert(v);
        }
        if (chosen.empty()) {
            // fallback: connect to a random existing node
            std::uniform_int_distribution<std::size_t> uni(0, u - 1);
            chosen.insert(uni(rng));
        }
        for (auto v : chosen) {
            double w = wr(rng);
            edges.emplace_back(u, v, w);
            if (undirected) edges.emplace_back(v, u, w);
            // update targets: increase degree of u and v
            std::size_t reps = (undirected ? 2 : 1);
            for (std::size_t t = 0; t < reps; ++t) targets.push_back(u);
            for (std::size_t t = 0; t < reps; ++t) targets.push_back(v);
        }
    }

    return build_csr(n, edges);
}

CSRGraph generate_rmat(int scale, int edgefactor, unsigned int seed,
                       double a, double b, double c, double wmin, double wmax, bool undirected) {
    std::size_t n = 1ull << scale;
    std::size_t M = static_cast<std::size_t>(edgefactor) * n;
    double d = 1.0 - (a + b + c);
    if (d < 0) d = 0.0;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> wr(wmin, wmax);

    std::vector<std::tuple<std::size_t,std::size_t,double>> edges;
    edges.reserve(M * (undirected ? 2 : 1));

    for (std::size_t e = 0; e < M; ++e) {
        std::size_t u = 0, v = 0;
        std::size_t bit = 1ull << (scale - 1);
        for (int s = 0; s < scale; ++s) {
            double r = urand(rng);
            if (r < a) {
                // (0,0)
            } else if (r < a + b) {
                v |= bit; // (0,1)
            } else if (r < a + b + c) {
                u |= bit; // (1,0)
            } else {
                u |= bit; v |= bit; // (1,1)
            }
            bit >>= 1;
        }
        if (u == v) { v = (v + 1) % n; }
        double w = wr(rng);
        edges.emplace_back(u, v, w);
        if (undirected) edges.emplace_back(v, u, w);
    }

    return build_csr(n, edges);
}
