#pragma once

#include <cstdint>

// Hash function selection enum
enum class HashFunctionType {
    STRONG,     // Strong hash (MurmurHash3-inspired) - good distribution
    WEAK,       // Weak hash (identity) - no scrambling, clusters preserved
    MODULO      // Modulo-only hash - minimal transformation
};

// Global hash function selector (default: STRONG)
extern HashFunctionType global_hash_function_type;

// ============= HASH FUNCTION IMPLEMENTATIONS =============

// Strong hash function: MurmurHash3-inspired 64-bit
// Ensures uniform distribution across hash space
inline uint64_t hash_strong(uint64_t key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key;
}

// Weak hash function: Identity (no scrambling)
// Preserves key patterns, causes clustering
inline uint64_t hash_weak(uint64_t key) {
    return key;
}

// Modulo-only hash function: Minimal transformation
// Uses only bit-level operations without multiplication
inline uint64_t hash_modulo(uint64_t key) {
    return (key ^ (key >> 16)) ^ (key >> 8);
}

// Dispatcher function: routes to selected hash function
inline uint64_t custom_hash(uint64_t key) {
    switch (global_hash_function_type) {
        case HashFunctionType::STRONG:
            return hash_strong(key);
        case HashFunctionType::WEAK:
            return hash_weak(key);
        case HashFunctionType::MODULO:
            return hash_modulo(key);
        default:
            return hash_strong(key);
    }
}