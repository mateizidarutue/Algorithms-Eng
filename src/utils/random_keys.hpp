#pragma once

#include <vector>
#include <random>
#include <cstdint>

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
