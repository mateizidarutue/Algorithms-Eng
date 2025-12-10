#pragma once

#include <vector>
#include <chrono>
#include <iostream>

#include "../hashing/robin_hood.hpp"

// Benchmarks our implementation of classic robin hood hashing
inline void benchmark_robinhood_custom(const std::vector<uint64_t>& keys, size_t capacity) {
    RobinHoodHashTable table(capacity);

    // Measure insertion time for classic manual implementation of Robin Hood
    auto start = std::chrono::high_resolution_clock::now();
    for (auto& key : keys) {
        table.insert(key);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "  Custom Robin Hood Insert Time: "
              << std::chrono::duration<double>(end - start).count() << " seconds\n";

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
    table.export_histograms_csv("rh");
}
