#pragma once

#include <vector>
#include <chrono>
#include <iostream>

#include "../hashing/robin_hood.hpp"
#include "../utils/random_keys.hpp"

// Benchmarks our implementation of classic robin hood hashing
inline void benchmark_robinhood_custom(const std::vector<uint64_t>& keys, size_t capacity, double load_factor = 0.0) {
    RobinHoodHashTable table(capacity);

    // Measure insertion time for classic manual implementation of Robin Hood
    auto start = std::chrono::high_resolution_clock::now();
    for (auto& key : keys) {
        table.insert(key);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "  Custom Robin Hood Insert Time: "
              << std::chrono::duration<double>(end - start).count() << " seconds\n";

    // Recompute final displacement histograms to reflect the completed table state
    table.recompute_final_distance_histograms();

    // HIT lookups
    start = std::chrono::high_resolution_clock::now();
    for (auto& key : keys) volatile bool f = table.lookup(key);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "  Custom Robin Hood Lookup HIT Time: "
              << std::chrono::duration<double>(end - start).count() << " s\n";

    // MISS lookups
    std::vector<uint64_t> miss_keys = generate_missing_keys(keys.size());

    start = std::chrono::high_resolution_clock::now();
    for (auto& key : miss_keys) volatile bool f = table.lookup(key);
    end = std::chrono::high_resolution_clock::now();

    std::cout << "  Custom Robin Hood Lookup MISS Time: "
              << std::chrono::duration<double>(end - start).count() << " s\n";

    table.print_probe_stats();
    
    // Create filename with load factor and hash function type
    std::string prefix = "rh";
    if (load_factor > 0.0) {
        prefix += "_lf_" + std::to_string(static_cast<int>(load_factor * 100));
    }
    table.export_histograms_csv(prefix);
}
