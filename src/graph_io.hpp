#pragma once
#include "graph.hpp"
#include <string>

// Load an edge list from a text file.
// Lines: optionally starting with '#' are comments and ignored.
// Format per edge: u v [w]
// - If has_weights is false, a default weight of 1.0 is used.
// - If one_based is true, vertex IDs are decremented by 1.
// - If undirected is true, edges are mirrored.
CSRGraph load_edge_list(const std::string& path,
                        bool has_weights = true,
                        bool one_based = true,
                        bool undirected = true);
