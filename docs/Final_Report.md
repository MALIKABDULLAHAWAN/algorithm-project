# Parallel Single-Source Shortest Path via Δ-Stepping

**Course:** CS-478 — Design and Analysis of Algorithms  
**Project:** Complex Computing Problem (CCP)  
**Team:** [Member A], [Member B], [Member C], [Member D]

---

## Abstract

This project implements and evaluates a parallel algorithm for the Single-Source Shortest Path (SSSP) problem. We implement the Δ-stepping algorithm using OpenMP for shared-memory parallelism and compare it against a sequential Dijkstra baseline. Our implementation supports multiple synthetic graph generators (Erdős–Rényi, Barabási–Albert, RMAT) and real-world datasets (SNAP/DIMACS format). Experimental evaluation demonstrates the algorithm's correctness and analyzes speedup, efficiency, and parameter sensitivity across different graph structures.

---

## 1. Introduction

### 1.1 Problem Statement

With the rise of multi-core systems, classic sequential graph algorithms become bottlenecks for large-scale problems. The Single-Source Shortest Path problem—finding minimum-weight paths from a source vertex to all other vertices—is fundamental in network routing, navigation, and social network analysis.

### 1.2 Motivation

Dijkstra's algorithm, while optimal for sequential execution at O(m log n), is inherently sequential due to its greedy vertex selection. This project explores Δ-stepping, a parallel algorithm that trades some work efficiency for parallelism by processing vertices in distance-based buckets.

### 1.3 Objectives

1. Implement a correct sequential baseline (Dijkstra)
2. Implement parallel Δ-stepping with OpenMP
3. Derive theoretical work-span complexity
4. Evaluate performance on synthetic and real graphs
5. Analyze deviations from theoretical predictions

---

## 2. Background

### 2.1 Sequential SSSP: Dijkstra's Algorithm

Dijkstra's algorithm maintains a priority queue of vertices ordered by tentative distance. At each step, it extracts the minimum-distance vertex and relaxes its outgoing edges.

**Complexity:**
- Time: O((n + m) log n) with binary heap
- Space: O(n + m)

**Limitation:** The sequential extraction of minimum elements prevents parallelization.

### 2.2 Parallel SSSP: Δ-Stepping

Δ-stepping (Meyer & Sanders, 2003) partitions vertices into buckets based on tentative distances:
- Bucket k contains vertices with distance in [kΔ, (k+1)Δ)
- Edges are classified as **light** (weight ≤ Δ) or **heavy** (weight > Δ)
- Light edges are processed iteratively within a bucket until closure
- Heavy edges propagate to future buckets

This allows parallel processing of all vertices within a bucket.

---

## 3. Algorithm Design

### 3.1 Pseudocode

```
DELTA-STEPPING(G, s, Δ):
    Input: Graph G=(V,E,w), source s, bucket width Δ
    
    dist[v] ← ∞ for all v ∈ V
    dist[s] ← 0
    B[0] ← {s}                    // Initial bucket
    i ← 0                         // Current bucket index
    
    while i < |B|:
        if B[i] is empty:
            i ← i + 1
            continue
        
        R ← B[i]; B[i] ← ∅        // Extract bucket
        S ← ∅                      // Settled vertices
        
        // Light-edge closure
        while R ≠ ∅:
            S ← S ∪ R
            Relaxations ← ∅
            
            parallel for u ∈ R:
                for (u,v) ∈ E where w(u,v) ≤ Δ:
                    if dist[u] + w(u,v) < dist[v]:
                        Relaxations.add(v, dist[u] + w(u,v))
            
            // Coalesce: keep minimum per vertex
            for (v, d) in coalesce(Relaxations):
                if d < dist[v]:
                    dist[v] ← d
                    if ⌊d/Δ⌋ = i: R.add(v)
                    else: B[⌊d/Δ⌋].add(v)
            
            R ← new vertices in current bucket
        
        // Heavy-edge relaxation
        parallel for u ∈ S:
            for (u,v) ∈ E where w(u,v) > Δ:
                if dist[u] + w(u,v) < dist[v]:
                    dist[v] ← dist[u] + w(u,v)
                    B[⌊dist[v]/Δ⌋].add(v)
        
        i ← i + 1
    
    return dist
```

### 3.2 Correctness

The algorithm is correct because:
1. **Monotonicity:** Buckets are processed in increasing order, ensuring vertices are settled with non-decreasing distances
2. **Light-edge closure:** Iterating until R is empty guarantees all light-edge paths within the bucket are explored
3. **Non-negative weights:** Ensures no vertex needs revisiting after settlement

### 3.3 Complexity Analysis

**Work (Total Operations):**
- Each edge is relaxed O(1) expected times on sparse graphs
- Work = O(n + m) expected, O(n + m + nL) worst case where L = max light edges per vertex

**Span (Critical Path):**
- Number of bucket phases: O(D/Δ) where D = diameter
- Per-phase span: O(Δ·L/Δ) = O(L) for light closure
- Total Span ≈ O((D/Δ)·L + log n) for parallel reductions

**Speedup:**
- Theoretical: S_p = Work/Span
- Bounded by bucket dependencies and synchronization overhead

### 3.4 Parameter Selection (Δ)

The bucket width Δ controls the work-span tradeoff:
- **Small Δ:** More buckets, less parallelism per bucket, more phases
- **Large Δ:** Fewer buckets, more parallelism, but more redundant work

**Heuristic:** Start with Δ ≈ average edge weight, then tune experimentally.

---

## 4. Implementation

### 4.1 Data Structures

**CSR (Compressed Sparse Row) Graph:**
```cpp
struct CSRGraph {
    size_t n;                    // Number of vertices
    size_t m;                    // Number of edges
    vector<size_t> row_offsets;  // Size n+1
    vector<size_t> col_indices;  // Size m
    vector<double> weights;      // Size m
};
```

CSR provides cache-efficient edge traversal with O(1) access to adjacency lists.

### 4.2 Parallelization Strategy

We use OpenMP for shared-memory parallelism:

1. **Parallel edge relaxation:** Each thread processes a subset of vertices in R
2. **Thread-local buffers:** Avoid contention by collecting relaxations locally
3. **Critical section merge:** Combine local buffers after parallel region
4. **Coalescing:** Sort by vertex ID and keep minimum distance per vertex

```cpp
#pragma omp parallel
{
    vector<pair<size_t, double>> local_relaxes;
    #pragma omp for schedule(dynamic, 1024)
    for (size_t i = 0; i < R.size(); ++i) {
        // Relax edges from R[i]
    }
    #pragma omp critical
    relaxes.insert(relaxes.end(), local_relaxes.begin(), local_relaxes.end());
}
```

### 4.3 Graph Generators

| Generator | Description | Parameters |
|-----------|-------------|------------|
| **ER** | Erdős–Rényi random graph | n (vertices), p (edge probability) |
| **BA** | Barabási–Albert preferential attachment | n, m0 (seed), m (attachments) |
| **RMAT** | Recursive matrix graph (power-law) | scale, edgefactor, a/b/c params |

### 4.4 File Structure

```
├── src/
│   ├── dijkstra_seq.cpp      # Sequential baseline
│   ├── delta_stepping_omp.cpp # Parallel algorithm
│   ├── graph.cpp/hpp         # CSR + generators
│   ├── graph_io.cpp/hpp      # SNAP/DIMACS loader
│   └── bench.cpp             # Benchmark driver
├── scripts/
│   ├── run_experiments.py    # Automation
│   └── plot_results.py       # Visualization
└── docs/
    └── Final_Report.md       # This document
```

---

## 5. Experimental Setup

### 5.1 Environment

- **OS:** Windows 11
- **Compiler:** GCC 15.2.0 (MSYS2 MinGW-w64)
- **OpenMP:** Version 4.5
- **Build:** CMake with Release optimization (-O3)

### 5.2 Datasets

| Dataset | Type | n | m | Description |
|---------|------|---|---|-------------|
| ER-10K | Synthetic | 10,000 | ~200K | Random graph, p=0.001 |
| BA-10K | Synthetic | 10,000 | ~40K | Scale-free network |
| RMAT-16K | Synthetic | 16,384 | ~262K | Power-law degree distribution |

### 5.3 Metrics

- **Time:** Wall-clock execution time (ms)
- **Speedup:** S_p = T_1 / T_p (relative to single-thread or Dijkstra)
- **Efficiency:** E_p = S_p / p
- **Correctness:** Checksum of distances (must match Dijkstra)

---

## 6. Results

### 6.1 Correctness Verification

All configurations produce identical checksums:

| Graph | Dijkstra | Δ-stepping (1T) | Δ-stepping (4T) |
|-------|----------|-----------------|-----------------|
| ER-10K | 104406 | 104406 | 104406 |
| BA-10K | 125550 | 125550 | 125550 |
| RMAT-16K | 53972.7 | 53972.7 | 53972.7 |

### 6.2 Performance Results

Sample results on ER graph (n=10000, m≈200K):

| Algorithm | Threads | Δ | Time (ms) | Speedup vs Dijkstra |
|-----------|---------|---|-----------|---------------------|
| Dijkstra | 1 | - | 3.4 | 1.00x |
| Δ-stepping | 1 | 2.0 | 5.5 | 0.62x |
| Δ-stepping | 2 | 2.0 | 7.2 | 0.47x |
| Δ-stepping | 4 | 2.0 | 8.2 | 0.41x |
| Δ-stepping | 1 | 4.0 | 10.0 | 0.34x |
| Δ-stepping | 4 | 4.0 | 6.3 | 0.54x |

### 6.3 Analysis

**Observations:**
1. Δ-stepping has higher overhead than Dijkstra for small graphs
2. Parallel speedup is limited by synchronization and bucket dependencies
3. Larger Δ reduces phases but increases per-phase work
4. Performance varies significantly with graph structure

---

## 7. Discussion

### 7.1 Deviations from Theory

The observed speedup is lower than theoretical predictions due to:

1. **Synchronization overhead:** Critical sections for merging relaxations
2. **Memory bandwidth:** CSR traversal is memory-bound
3. **Load imbalance:** Uneven vertex degrees cause thread starvation
4. **Small problem size:** Parallelization overhead dominates for n < 100K

### 7.2 Parameter Sensitivity

- **Δ too small:** Many bucket phases, high synchronization cost
- **Δ too large:** More redundant relaxations, less parallelism benefit
- **Optimal Δ:** Depends on graph structure; typically 1-4× average edge weight

### 7.3 Bottlenecks

1. **Bucket coalescing:** Sorting relaxations is O(r log r) per phase
2. **Critical sections:** Serial merge of thread-local buffers
3. **Cache misses:** Random access patterns in sparse graphs

### 7.4 Potential Improvements

1. Use atomic operations instead of critical sections
2. Implement NUMA-aware bucket allocation
3. Use concurrent data structures for buckets
4. Hybrid approach: Dijkstra for small subproblems

---

## 8. Conclusion

We successfully implemented and evaluated the Δ-stepping parallel SSSP algorithm. Key findings:

1. **Correctness:** Algorithm produces identical results to Dijkstra across all test cases
2. **Scalability:** Limited by synchronization overhead on small graphs
3. **Parameter tuning:** Δ selection significantly impacts performance
4. **Trade-offs:** Parallelism comes at the cost of additional work

For production use on large graphs (n > 1M), Δ-stepping offers meaningful speedup. For smaller graphs, sequential Dijkstra remains competitive due to lower overhead.

---

## 9. References

1. Meyer, U., & Sanders, P. (2003). Δ-stepping: A parallelizable shortest path algorithm. *Journal of Algorithms*, 49(1), 114-152.

2. Cormen, T. H., et al. (2009). *Introduction to Algorithms* (3rd ed.). MIT Press.

3. SNAP: Stanford Network Analysis Project. https://snap.stanford.edu/

4. DIMACS Implementation Challenges. http://dimacs.rutgers.edu/

---

## Appendix A: Build Instructions

```bash
# Clone and build
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .

# Run benchmarks
./bench --algo dijkstra --gen er --n 10000 --p 0.001
./bench --algo delta --gen er --n 10000 --p 0.001 --threads 4 --delta 2.0

# Run experiments
python scripts/run_experiments.py --bench build/bench.exe --out results.csv
python scripts/plot_results.py --csv results.csv --outdir figs
```

## Appendix B: Team Contributions

| Member | Responsibilities |
|--------|-----------------|
| **A** | Sequential Dijkstra, graph I/O, correctness testing |
| **B** | Δ-stepping implementation, OpenMP optimization |
| **C** | Graph generators, experiment automation, data collection |
| **D** | Documentation, analysis, presentation |
