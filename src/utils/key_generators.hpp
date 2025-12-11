#pragma once

#include <vector>
#include <cstdint>
#include <random>
#include <algorithm>

class KeyGenerator {
public:
    virtual ~KeyGenerator() = default;
    virtual std::vector<uint64_t> generate(size_t count, uint64_t seed = 42) = 0;
};

class RandomKeyGenerator : public KeyGenerator {
public:
    std::vector<uint64_t> generate(size_t count, uint64_t seed = 42) override {
        std::vector<uint64_t> keys(count);
        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<uint64_t> dist;
        
        for (size_t i = 0; i < count; ++i) {
            keys[i] = dist(rng);
        }
        return keys;
    }
};

class ClusteredKeyGenerator : public KeyGenerator {
public:
    std::vector<uint64_t> generate(size_t count, uint64_t seed = 42) override {
        std::vector<uint64_t> keys(count);
        std::mt19937_64 rng(seed);
        
        // Generate clustered keys: many keys in tight ranges
        // This creates clustering in hash space when mapped to table indices
        size_t num_clusters = (count / 100) + 1;  // ~100 keys per cluster
        size_t keys_per_cluster = count / num_clusters;
        
        size_t idx = 0;
        for (size_t cluster = 0; cluster < num_clusters && idx < count; ++cluster) {
            // Each cluster has a base value
            uint64_t base = rng() & 0xFFFFFFFF00000000ULL;  // High 32 bits define cluster
            
            // Generate small offsets within cluster
            std::uniform_int_distribution<uint64_t> offset_dist(0, 0xFFFFULL);  // Low 16 bits
            
            for (size_t i = 0; i < keys_per_cluster && idx < count; ++i) {
                keys[idx++] = base | offset_dist(rng);
            }
        }
        
        return keys;
    }
};

class AdversarialKeyGenerator : public KeyGenerator {
public:
    std::vector<uint64_t> generate(size_t count, uint64_t seed = 42) override {
        std::vector<uint64_t> keys(count);
        
        // Generate adversarial keys: sequential integers that collide under weak hash
        // Many keys will hash to nearby slots, causing clustering
        for (size_t i = 0; i < count; ++i) {
            keys[i] = static_cast<uint64_t>(i);
        }
        
        return keys;
    }
};

class BitReversedKeyGenerator : public KeyGenerator {
public:
    std::vector<uint64_t> generate(size_t count, uint64_t seed = 42) override {
        std::vector<uint64_t> keys(count);
        
        // Generate sequential keys and bit-reverse them
        // This creates patterns that might collide with certain hash functions
        for (size_t i = 0; i < count; ++i) {
            keys[i] = reverse_bits(static_cast<uint64_t>(i));
        }
        
        return keys;
    }

private:
    static uint64_t reverse_bits(uint64_t x) {
        x = ((x & 0xAAAAAAAAAAAAAAAAULL) >> 1) | ((x & 0x5555555555555555ULL) << 1);
        x = ((x & 0xCCCCCCCCCCCCCCCCULL) >> 2) | ((x & 0x3333333333333333ULL) << 2);
        x = ((x & 0xF0F0F0F0F0F0F0F0ULL) >> 4) | ((x & 0x0F0F0F0F0F0F0F0FULL) << 4);
        x = ((x & 0xFF00FF00FF00FF00ULL) >> 8) | ((x & 0x00FF00FF00FF00FFULL) << 8);
        x = ((x & 0xFFFF0000FFFF0000ULL) >> 16) | ((x & 0x0000FFFF0000FFFFULL) << 16);
        x = (x >> 32) | (x << 32);
        return x;
    }
};

// Generate missing keys (keys not in the inserted set)
inline std::vector<uint64_t> generate_missing_keys(const std::vector<uint64_t>& inserted_keys, size_t count, uint64_t seed = 43) {
    std::vector<uint64_t> missing;
    missing.reserve(count);
    
    std::set<uint64_t> inserted_set(inserted_keys.begin(), inserted_keys.end());
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<uint64_t> dist;
    
    while (missing.size() < count) {
        uint64_t key = dist(rng);
        if (inserted_set.find(key) == inserted_set.end()) {
            missing.push_back(key);
            inserted_set.insert(key);  // Track to avoid duplicates
        }
    }
    
    return missing;
}
