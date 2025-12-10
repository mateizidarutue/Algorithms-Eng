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
    const size_t table_size = 1'000'000;

    // CONFIGURABLE: Number of keys to insert based on load factor
    const size_t num_keys_to_insert = static_cast<size_t>(table_size * load_factor);

    // Generate random keys
    std::vector<uint64_t> keys_random = generate_random_keys(num_keys_to_insert);

    std::cout << "\033[31mBENCHMARKING WIH RANDOM KEYS:\033[0m\n";
    //benchmark_robin_hood_library(keys); // tsl::robin_map
    benchmark_linear_probing(keys_random, table_size);
    benchmark_robinhood_custom(keys_random, table_size);

    // Generate clustered keys
    std::vector<uint64_t> keys_clustered = generate_clustered_keys(num_keys_to_insert);

    // Benchmark with clustered keys
    std::cout << "\033[31mBENCHMARKING WITH CLUSTERED KEYS:\033[0m\n";
    benchmark_linear_probing(keys_clustered, table_size);
    benchmark_robinhood_custom(keys_clustered, table_size);

    return 0;
}
