#pragma once

#include <vector>
#include <chrono>
#include <iostream>

#include "../hashing/linear_probing.hpp"

// Encapsulates linear probing insertion and lookup timing
inline void benchmark_linear_probing(const std::vector<uint64_t>& keys, size_t table_capacity) {
    // Linear Probing Hash Table
    LinearProbingHashTable table(table_capacity);

    // Measure insertion time for Linear Probing
    auto start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) {
        table.insert(key);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Linear Probing Insert Time: "
              << std::chrono::duration<double>(end - start).count() << " seconds\n";

    // Measure lookup time for Linear Probing
    start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) {
        volatile bool found = table.lookup(key);
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Linear Probing Lookup Time: "
              << std::chrono::duration<double>(end - start).count() << " seconds\n";

    table.print_probe_stats();
    table.export_histograms_csv("lp");
}
