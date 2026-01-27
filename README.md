# Parallel SSSP via Δ-Stepping

CS-478 Complex Computing Problem: Parallel Single-Source Shortest Path implementation using OpenMP.

## Quick Start

```bash
# Build
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .

# Run
./bench --algo dijkstra --gen er --n 10000 --p 0.001
./bench --algo delta --gen er --n 10000 --p 0.001 --threads 4 --delta 2.0
```

## Features

- Sequential Dijkstra baseline (O(m log n))
- Parallel Δ-stepping with OpenMP
- Graph generators: ER, BA, RMAT
- SNAP/DIMACS edge-list loader
- Automated experiments and plotting

## Documentation

See [docs/Final_Report.md](docs/Final_Report.md) for complete analysis including:
- Algorithm design and pseudocode
- Work-span complexity analysis
- Experimental results
- Performance discussion

## Project Structure

```
src/           - C++ implementation
scripts/       - Python experiment/plotting scripts
docs/          - Final report
build/         - Compiled binaries
```

## Usage

```bash
# Experiments
python scripts/run_experiments.py --bench build/bench.exe --out results.csv

# Plots
python scripts/plot_results.py --csv results.csv --outdir figs
```
