#pragma once
#include <vector>
#include <cstddef>
#include <random>
#include <utility>

struct CSRGraph {
    using index_t = std::size_t;
    using weight_t = double;

    index_t n = 0;              // number of vertices
    index_t m = 0;              // number of edges (directed)
    std::vector<index_t> row_offsets; // size n+1
    std::vector<index_t> col_indices; // size m
    std::vector<weight_t> weights;    // size m

    bool empty() const { return n == 0; }
};

// Build CSR from edge list (0-based u,v, directed). Will coalesce into CSR without sorting.
CSRGraph build_csr(std::size_t n, const std::vector<std::tuple<std::size_t,std::size_t,double>>& edges);

// Generate an Erdős–Rényi G(n, p) directed graph with weights in [wmin, wmax].
CSRGraph generate_er(std::size_t n, double p, unsigned int seed,
                     double wmin = 1.0, double wmax = 10.0,
                     bool undirected = true);

// Generate a Barabási–Albert preferential attachment graph.
// Start with m0 fully connected seed nodes, then add nodes each with m_attach edges.
CSRGraph generate_ba(std::size_t n, std::size_t m0, std::size_t m_attach, unsigned int seed,
                     double wmin = 1.0, double wmax = 10.0,
                     bool undirected = true);

// Generate an RMAT graph with 2^scale nodes and approximately edgefactor * 2^scale edges.
// RMAT parameters a,b,c (with d=1-a-b-c). Typical: a=0.57, b=0.19, c=0.19.
CSRGraph generate_rmat(int scale, int edgefactor, unsigned int seed,
                       double a = 0.57, double b = 0.19, double c = 0.19,
                       double wmin = 1.0, double wmax = 10.0,
                       bool undirected = true);
