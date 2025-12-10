#pragma once

#include <vector>
#include <random>
#include <cstdint>

#include "hashing/hash_function.hpp"

inline std::vector<uint64_t> generate_random_keys(size_t n, uint64_t seed = 42) {
    std::vector<uint64_t> keys(n);

    // CONFIGURABLE: Randomness source, seed, distribution, duplicate-key frequency
    // Impacts clustering and collision patterns
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<uint64_t> dist;

    for (size_t i = 0; i < n; i++)
        keys[i] = dist(rng);

    return keys;
}

// Clustered keys: simulate real-world locality
// Effect dimished by powerful custom_hash function inside hashing/hash_function.hpp
inline std::vector<uint64_t> generate_clustered_keys(size_t n) { 
    std::vector<uint64_t> keys(n); 
    std::mt19937_64 rng(1337); 

    const size_t clusters = 10; 
    const size_t cluster_size = n / clusters; 

    for (size_t c = 0; c < clusters; ++c) { 
        uint64_t base = 1'000'000ULL * (c + 1); 
        std::uniform_int_distribution<uint64_t> offset(0, 5000); 

        for (size_t i = 0; i < cluster_size; ++i) { 
            keys[c * cluster_size + i] = base + offset(rng); 
        } 
    } 
    return keys;
}

// Generates multiple clusters of keys that hash into specific regions of the table
// Takes a lot of time to run for large input
inline std::vector<uint64_t> generate_multi_cluster_keys(
        size_t n, size_t capacity, size_t clusters = 8, size_t cluster_width = 50)
{
    std::vector<uint64_t> keys;
    keys.reserve(n);

    size_t keys_per_cluster = n / clusters;

    uint64_t candidate = 0;
    for (size_t c = 0; c < clusters; c++) {
        size_t bucket_min = c * cluster_width;
        size_t bucket_max = bucket_min + cluster_width;

        size_t count = 0;
        while (count < keys_per_cluster) {
            size_t bucket = custom_hash(candidate) % capacity;
            if (bucket >= bucket_min && bucket < bucket_max) {
                keys.push_back(candidate);
                count++;
            }
            candidate++;
        }
    }
    return keys;
}


// Keys guaranteed NOT to exist in the table = used for MISS lookups
inline std::vector<uint64_t> generate_missing_keys(size_t n) {
    std::vector<uint64_t> miss = generate_random_keys(n, 99999);

    for (auto& k : miss)
        k += 5'000'000'000ULL;  // ensure no overlap with inserted keys

    return miss;
}