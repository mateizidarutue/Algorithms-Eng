#pragma once

#include <vector>
#include <chrono>
#include <iostream>

#include "../hashing/linear_probing.hpp"

// Encapsulates linear probing insertion and lookup timing
inline void benchmark_linear_probing(const std::vector<uint64_t>& keys, size_t table_capacity, double load_factor = 0.0) {
    // Linear Probing Hash Table
    LinearProbingHashTable table(table_capacity);

    // Measure insertion time for Linear Probing
    auto start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) {
        table.insert(key);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "  Linear Probing Insert Time: "
              << std::chrono::duration<double>(end - start).count() << " seconds\n";

    // HIT lookups = lookup keys we actually inserted
    start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) volatile bool f = table.lookup(key);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "  Linear Probing Lookup HIT Time: "
              << std::chrono::duration<double>(end - start).count() << " s\n";

    // MISS lookups = lookup keys that are not present
    std::vector<uint64_t> miss_keys = generate_missing_keys(keys.size());

    start = std::chrono::high_resolution_clock::now();
    for (const auto& key : miss_keys) volatile bool f = table.lookup(key);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "  Linear Probing Lookup MISS Time: "
              << std::chrono::duration<double>(end - start).count() << " s\n";

    table.print_probe_stats();
    
    // Create filename with load factor and hash function type
    std::string prefix = "lp";
    if (load_factor > 0.0) {
        prefix += "_lf_" + std::to_string(static_cast<int>(load_factor * 100));
    }
    table.export_histograms_csv(prefix);
}
