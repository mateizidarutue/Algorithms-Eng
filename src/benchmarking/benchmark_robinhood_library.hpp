// #pragma once

// #include <vector>
// #include <tsl/robin_map.h>
// #include <chrono>
// #include <iostream>

// #include "../hashing/hash_function.hpp"

// // Encapsulates library implementation of Robin Hood hashing: insertion and lookup timing
// inline void benchmark_robin_hood_library(const std::vector<uint64_t>& keys) {
//     // Robin Hood Hashing
//     tsl::robin_map<uint64_t, bool, decltype(&custom_hash)> robin_map(0, custom_hash);

//     // Measure insertion time for Robin Hood Hashing
//     auto start = std::chrono::high_resolution_clock::now();
//     for (const auto& key : keys) {
//         robin_map[key] = true;
//     }
//     auto end = std::chrono::high_resolution_clock::now();
//     std::cout << "Robin Hood (library) Insert Time: "
//               << std::chrono::duration<double>(end - start).count() << " seconds\n";

//     // Measure lookup time for Robin Hood Hashing
//     start = std::chrono::high_resolution_clock::now();
//     for (const auto& key : keys) {
//         volatile bool found = robin_map.find(key) != robin_map.end();
//     }
//     end = std::chrono::high_resolution_clock::now();
//     std::cout << "Robin Hood (library) Lookup Time: "
//               << std::chrono::duration<double>(end - start).count() << " seconds\n";
// }

