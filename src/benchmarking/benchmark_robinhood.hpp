#pragma once

#include <vector>
#include <chrono>
#include <iostream>
#include <string>
#include <sstream>

#include "../hashing/robin_hood.hpp"
#include "../utils/random_keys.hpp"

// Benchmarks our implementation of classic robin hood hashing
inline void benchmark_robinhood_custom(const std::vector<uint64_t>& keys,
                                       size_t capacity,
                                       double load_factor,
                                       const std::string& dataset_label) {
    RobinHoodHashTable table(capacity);

    // Measure insertion time for classic manual implementation of Robin Hood
    auto start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) {
        table.insert(key);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "  Custom Robin Hood Insert Time: "
              << std::chrono::duration<double>(end - start).count()
              << " seconds\n";

    // HIT lookups
    start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) {
        volatile bool f = table.lookup(key);
        (void)f;
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "  Custom Robin Hood Lookup HIT Time: "
              << std::chrono::duration<double>(end - start).count()
              << " s\n";

    // MISS lookups
    std::vector<uint64_t> miss_keys = generate_missing_keys(keys.size());

    start = std::chrono::high_resolution_clock::now();
    for (const auto& key : miss_keys) {
        volatile bool f = table.lookup(key);
        (void)f;
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "  Custom Robin Hood Lookup MISS Time: "
              << std::chrono::duration<double>(end - start).count()
              << " s\n";

    // Print stats
    table.print_probe_stats();

    // Export histograms with prefix that encodes LF + dataset type
    std::ostringstream oss;
    int lf_percent = static_cast<int>(load_factor * 100.0 + 0.5);
    oss << "rh_lf_" << lf_percent << "_" << dataset_label;

    table.export_histograms_csv(oss.str());
}
