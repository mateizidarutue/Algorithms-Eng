#include <iostream>
#include "utils/random_keys.hpp"
#include "benchmarking/benchmark_robinhood_library.hpp"
#include "benchmarking/benchmark_linearprobing.hpp"
#include "benchmarking/benchmark_robinhood.hpp"

int main() {
    /* CONFIGURABLE:
        Load factor: number of elements inserted / table capacity
        Higher load -> more collisions -> more probes
        Lower load factor -> faster lookups, fewer collisions
    */
    const double load_factor = 0.7;

    // CONFIGURABLE: Table size: number of slots in the hash table for testing different behaviors
    const size_t num_elements = 1'000'000;

    // CONFIGURABLE: Number of keys to insert based on load factor
    const size_t num_keys_to_insert = static_cast<size_t>(num_elements * load_factor);

    // Generate random keys
    std::vector<uint64_t> keys = generate_random_keys(num_keys_to_insert);

    //benchmark_robin_hood_library(keys); // tsl::robin_map
    benchmark_linear_probing(keys, num_elements);
    benchmark_robinhood_custom(keys, num_elements);

    return 0;
}
