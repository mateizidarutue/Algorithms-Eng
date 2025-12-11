#include <iostream>
#include <algorithm>
#include "utils/random_keys.hpp"
#include "benchmarking/benchmark_robinhood_library.hpp"
#include "benchmarking/benchmark_linearprobing.hpp"
#include "benchmarking/benchmark_robinhood.hpp"
#include "hashing/hash_function.hpp"

// Define the global hash function type
HashFunctionType global_hash_function_type = HashFunctionType::WEAK;

int main() {
    // CONFIGURABLE: Table size: number of slots in the hash table for testing different behaviors
    const size_t table_size = 1'000'000;

    /* CONFIGURABLE:
        Load factors to test: number of elements inserted / table capacity
        Higher load -> more collisions -> more probes
        Lower load factor -> faster lookups, fewer collisions
    */
    std::vector<double> load_factors = {
        0.25, 0.40, 0.50, 0.70, 0.80, 0.90,
        0.95, 0.96, 0.97, 0.98, 0.99
    };

    const double max_lf = *std::max_element(load_factors.begin(), load_factors.end());
    const size_t max_keys = static_cast<size_t>(table_size * max_lf);

    // Shared supersets so runs are comparable across load factors
    std::vector<uint64_t> all_keys_random = generate_random_keys(max_keys);
    std::vector<uint64_t> all_keys_clustered = generate_clustered_keys(max_keys);

    for (double load_factor : load_factors) {
        std::cout << "\n\033[1;33m========== Testing Load Factor: " << load_factor << " ==========\033[0m" << std::endl;
        
        // CONFIGURABLE: Number of keys to insert based on load factor
        const size_t num_keys_to_insert = static_cast<size_t>(table_size * load_factor);

        // Slice from the shared supersets for reproducibility across LFs
        std::vector<uint64_t> keys_random(all_keys_random.begin(),
                                          all_keys_random.begin() + num_keys_to_insert);

        std::cout << "\033[31mBENCHMARKING WITH RANDOM KEYS:\033[0m\n";
        //benchmark_robin_hood_library(keys); // tsl::robin_map
        benchmark_linear_probing(keys_random, table_size, load_factor);
        benchmark_robinhood_custom(keys_random, table_size, load_factor);

        // Generate clustered keys from the same superset
        //std::vector<uint64_t> keys_clustered(all_keys_clustered.begin(),
        //                                     all_keys_clustered.begin() + num_keys_to_insert);

        // Benchmark with clustered keys
        //std::cout << "\033[31mBENCHMARKING WITH CLUSTERED KEYS:\033[0m\n";
        //benchmark_linear_probing(keys_clustered, table_size, load_factor);
        //benchmark_robinhood_custom(keys_clustered, table_size, load_factor);
    }

    return 0;
}
