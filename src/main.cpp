#include <iostream>
#include <vector>

#include "utils/random_keys.hpp"
#include "benchmarking/benchmark_robinhood_library.hpp"
#include "benchmarking/benchmark_linearprobing.hpp"
#include "benchmarking/benchmark_robinhood.hpp"

int main() {
    const size_t table_size = 1'000'000;

    // Load factors we want to test
    const std::vector<double> load_factors = {
        0.25, 0.40, 0.50, 0.70, 0.80, 0.90, 0.95, 0.97, 0.99
    };

    for (double lf : load_factors) {
        size_t num_keys = static_cast<size_t>(lf * table_size);

        // ---------- RANDOM KEYS ----------
        auto keys_random = generate_random_keys(num_keys);

        std::cout << "\n=== LF = " << lf << " (RANDOM) ===\n";
        benchmark_linear_probing(keys_random, table_size, lf, "random");
        benchmark_robinhood_custom(keys_random, table_size, lf, "random");

        // ---------- MULTI-CLUSTER (HASH-CLUSTERED) KEYS ----------
        auto keys_multi = generate_multi_cluster_keys(num_keys, table_size);

        std::cout << "\n=== LF = " << lf << " (MULTI-CLUSTER) ===\n";
        benchmark_linear_probing(keys_multi, table_size, lf, "multi");
        benchmark_robinhood_custom(keys_multi, table_size, lf, "multi");
    }

    return 0;
}
