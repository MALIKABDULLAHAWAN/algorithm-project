#include "graph.hpp"
#include "graph_io.hpp"
#include <chrono>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#ifdef HAS_OPENMP
#include <omp.h>
#endif

// forward declarations
std::vector<CSRGraph::weight_t> dijkstra_seq(const CSRGraph& g, std::size_t source);
std::vector<CSRGraph::weight_t> delta_stepping_omp(const CSRGraph& g, std::size_t source, double delta, int threads);

struct Args {
    std::string algo = "dijkstra"; // or "delta"
    std::string gen = "er"; // er|ba|rmat when no input file
    std::size_t n = 10000;
    double p = 0.0005; // ER probability
    unsigned seed = 42;
    std::size_t source = 0;
    int threads = 1;
    double delta = 1.0;
    std::string input = ""; // path to edge list
    bool has_weights = true;
    bool one_based = true;
    bool undirected = true;
    std::string csv = ""; // output csv path
    // BA params
    std::size_t m0 = 4;
    std::size_t m_attach = 2;
    // RMAT params
    int scale = 16;
    int edgefactor = 16;
    double a = 0.57, b = 0.19, c = 0.19;
};

Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        auto next = [&](const char* name){ if (i + 1 >= argc) { std::cerr << "Missing value for " << name << "\n"; std::exit(1);} return std::string(argv[++i]); };
        if (k == "--algo") a.algo = next("--algo");
        else if (k == "--gen") a.gen = next("--gen");
        else if (k == "--n") a.n = std::stoull(next("--n"));
        else if (k == "--p") a.p = std::stod(next("--p"));
        else if (k == "--seed") a.seed = static_cast<unsigned>(std::stoul(next("--seed")));
        else if (k == "--source") a.source = std::stoull(next("--source"));
        else if (k == "--threads") a.threads = std::stoi(next("--threads"));
        else if (k == "--delta") a.delta = std::stod(next("--delta"));
        else if (k == "--input") a.input = next("--input");
        else if (k == "--weights") a.has_weights = (std::stoi(next("--weights")) != 0);
        else if (k == "--one_based") a.one_based = (std::stoi(next("--one_based")) != 0);
        else if (k == "--undirected") a.undirected = (std::stoi(next("--undirected")) != 0);
        else if (k == "--csv") a.csv = next("--csv");
        else if (k == "--m0") a.m0 = std::stoull(next("--m0"));
        else if (k == "--m_attach") a.m_attach = std::stoull(next("--m_attach"));
        else if (k == "--scale") a.scale = std::stoi(next("--scale"));
        else if (k == "--edgefactor") a.edgefactor = std::stoi(next("--edgefactor"));
        else if (k == "--a") a.a = std::stod(next("--a"));
        else if (k == "--b") a.b = std::stod(next("--b"));
        else if (k == "--c") a.c = std::stod(next("--c"));
        else if (k == "--help" || k == "-h") {
            std::cout << "Usage: bench [--algo dijkstra|delta] [--gen er|ba|rmat] [--n N] [--p P] [--seed S] [--source S0] [--threads T] [--delta D]\\n"
                         "            [--input FILE] [--weights 0|1] [--one_based 0|1] [--undirected 0|1] [--csv FILE]\\n"
                         "            [--m0 M0] [--m_attach M] [--scale S] [--edgefactor EF] [--a A] [--b B] [--c C]\n";
            std::exit(0);
        }
    }
    return a;
}

int main(int argc, char** argv) {
    auto args = parse_args(argc, argv);

#ifdef HAS_OPENMP
    if (args.threads > 0) {
        omp_set_num_threads(args.threads);
    }
#endif

    CSRGraph g;
    if (!args.input.empty()) {
        std::cout << "Loading graph from: " << args.input << " (weights=" << (args.has_weights?1:0)
                  << ", one_based=" << (args.one_based?1:0) << ", undirected=" << (args.undirected?1:0) << ")\n";
        try {
            g = load_edge_list(args.input, args.has_weights, args.one_based, args.undirected);
        } catch (const std::exception& e) {
            std::cerr << e.what() << "\n";
            return 1;
        }
    } else {
        if (args.gen == "er") {
            std::cout << "Generating ER graph: n=" << args.n << ", p=" << args.p << ", seed=" << args.seed << "\n";
            g = generate_er(args.n, args.p, args.seed, 1.0, 10.0, args.undirected);
        } else if (args.gen == "ba") {
            std::cout << "Generating BA graph: n=" << args.n << ", m0=" << args.m0 << ", m_attach=" << args.m_attach << ", seed=" << args.seed << "\n";
            g = generate_ba(args.n, args.m0, args.m_attach, args.seed, 1.0, 10.0, args.undirected);
        } else if (args.gen == "rmat") {
            std::cout << "Generating RMAT graph: scale=" << args.scale << ", edgefactor=" << args.edgefactor << ", seed=" << args.seed
                      << ", a=" << args.a << ", b=" << args.b << ", c=" << args.c << "\n";
            g = generate_rmat(args.scale, args.edgefactor, args.seed, args.a, args.b, args.c, 1.0, 10.0, args.undirected);
        } else {
            std::cerr << "Unknown generator: " << args.gen << "\n";
            return 1;
        }
    }
    if (args.source >= g.n) args.source = 0;

    std::cout << "Graph stats: n=" << g.n << ", m=" << g.m << " (directed edges)\n";

    auto t0 = std::chrono::high_resolution_clock::now();
    std::vector<CSRGraph::weight_t> dist;
    if (args.algo == "dijkstra") {
        dist = dijkstra_seq(g, args.source);
    } else if (args.algo == "delta") {
        dist = delta_stepping_omp(g, args.source, args.delta, args.threads);
    } else {
        std::cerr << "Unknown algo: " << args.algo << "\n";
        return 1;
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Basic checksum to avoid optimizing away
    double sum = 0.0; int cnt = 0;
    for (auto d : dist) { if (d < 1e300) { sum += d; ++cnt; } }

    std::cout << "Algo=" << args.algo
              << ", threads=" << args.threads
              << ", delta=" << args.delta
              << ", time_ms=" << ms
              << ", reachable=" << cnt
              << ", checksum=" << sum << "\n";

    if (!args.csv.empty()) {
        bool exists = std::filesystem::exists(args.csv);
        std::ofstream fout(args.csv, std::ios::app);
        if (!exists) {
            fout << "algo,threads,delta,n,m,time_ms,reachable,checksum,input,weights,one_based,undirected,gen,p,m0,m_attach,scale,edgefactor,a,b,c\n";
        }
        std::string input_str = args.input.empty() ? args.gen : args.input;
        fout << args.algo << "," << args.threads << "," << args.delta << ","
             << g.n << "," << g.m << "," << ms << "," << cnt << "," << sum << ","
             << "\"" << input_str << "\"" << ","
             << (args.has_weights?1:0) << "," << (args.one_based?1:0) << "," << (args.undirected?1:0)
             << "," << args.gen << "," << args.p << "," << args.m0 << "," << args.m_attach
             << "," << args.scale << "," << args.edgefactor << "," << args.a << "," << args.b << "," << args.c
             << "\n";
    }

    return 0;
}
