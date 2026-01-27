#include "graph_io.hpp"
#include <fstream>
#include <sstream>
#include <tuple>
#include <vector>
#include <cctype>

static inline bool starts_with_hash(const std::string& s) {
    for (char c : s) {
        if (std::isspace(static_cast<unsigned char>(c))) continue;
        return c == '#';
    }
    return false;
}

CSRGraph load_edge_list(const std::string& path, bool has_weights, bool one_based, bool undirected) {
    std::ifstream fin(path);
    if (!fin) {
        throw std::runtime_error("Cannot open file: " + path);
    }
    std::vector<std::tuple<std::size_t,std::size_t,double>> edges;
    edges.reserve(1024);
    std::size_t max_id = 0;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || starts_with_hash(line)) continue;
        std::istringstream iss(line);
        long long u=-1, v=-1; double w=1.0;
        if (has_weights) {
            if (!(iss >> u >> v >> w)) continue;
        } else {
            if (!(iss >> u >> v)) continue;
        }
        if (one_based) { u -= 1; v -= 1; }
        if (u < 0 || v < 0) continue;
        auto uu = static_cast<std::size_t>(u);
        auto vv = static_cast<std::size_t>(v);
        edges.emplace_back(uu, vv, w);
        if (undirected) edges.emplace_back(vv, uu, w);
        max_id = std::max(max_id, std::max<std::size_t>(uu, vv));
    }
    return build_csr(max_id + 1, edges);
}
