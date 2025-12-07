#include <iostream>
#include <vector>
#include <tsl/robin_map.h>
#include <chrono>
#include <random>

// Hash function for reproducibility ensuring uniform distribution of hash outputs
uint64_t custom_hash(uint64_t key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key;
}


class LinearProbingHashTable
{
public:
    explicit LinearProbingHashTable(size_t size)
        : table(size,  {false, 0}), capacity(size) {}

    void insert(uint64_t key) {
        size_t index = custom_hash(key) % capacity;

        while (table[index].occupied) {
            if (table[index].key == key) {
                return; // Key already exists
            }

            index = (index + 1) % capacity;
        }

        table[index] = {true, key};
    }

    // Checks if a given key exists in the hash table.
    // This function uses linear probing to search for the specified key
    // in the hash table. It starts at the index determined by the hash
    // of the key and iterates through the table until it finds the key
    // or encounters an unoccupied slot.
    bool lookup(uint64_t key) const {
        size_t index = custom_hash(key) % capacity;

        // Linear probing to find the key 
        while (table[index].occupied) {
            if (table[index].key == key) {
                return true; 
            }

            index = (index + 1) % capacity;
        }

        return false;
    }

private:
    // Represents a slot in the hash table, indicating whether it is occupied and storing the key if it is.
    struct Slot
    {
        bool occupied;
        uint64_t key;
    };
    std::vector<Slot> table; // The hash table storing slots: pair of (occupied flag, key)
    size_t capacity;
};

// Encapsulates Robin Hood hashing insertion and lookup timing
void benchmark_robin_hood(const std::vector<uint64_t>& keys) {
    // Robin Hood Hashing
    tsl::robin_map<uint64_t, bool, decltype(&custom_hash)> robin_map(0, custom_hash);
    
    // Measure insertion time for Robin Hood Hashing
    auto start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) {
        robin_map[key] = true;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> robin_insert_duration = end - start;
    std::cout << "Robin Hood Insertion Time: " << robin_insert_duration.count() << " seconds\n";

    // Measure lookup time for Robin Hood Hashing
    start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) {
        volatile bool found = robin_map.find(key) != robin_map.end();
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> robin_lookup_duration = end - start;
    std::cout << "Robin Hood Lookup Time: " << robin_lookup_duration.count() << " seconds\n";
}

// Encapsulates linear probing insertion and lookup timing
void benchmark_linear_probing(const std::vector<uint64_t>& keys, size_t table_capacity) {
    // Linear Probing Hash Table
    LinearProbingHashTable linear_probing_table(table_capacity);

    // Measure insertion time for Linear Probing
    auto start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) {
        linear_probing_table.insert(key);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> linear_insert_duration = end - start;
    std::cout << "Linear Probing Insertion Time: " << linear_insert_duration.count() << " seconds\n";

    // Measure lookup time for Linear Probing
    start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) {
        volatile bool found = linear_probing_table.lookup(key);
    }
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> linear_lookup_duration = end - start;
    std::cout << "Linear Probing Lookup Time: " << linear_lookup_duration.count() << " seconds\n";
}

int main() {
    /* Load factor: number of elements inserted / table capacity
       Higher load results more collisions and probing. */
    const double load_factor = 0.7; 
    const size_t num_elements = 1'000'000;
    const size_t num_keys_to_insert = static_cast<size_t>(num_elements * load_factor);

    // Generate random keys
    std::vector<uint64_t> keys(num_keys_to_insert);
    std::mt19937_64 rng(42); // Fixed seed for reproducibility
    std::uniform_int_distribution<uint64_t> dist;

    // Generates a set of random keys to be inserted into the hash table.
    // The generated keys will be used for testing the tame it takes to insert and lookup keys
    for (size_t i = 0; i < num_keys_to_insert; ++i) {
        keys[i] = dist(rng);
    }

    benchmark_robin_hood(keys);
    benchmark_linear_probing(keys, num_elements);

    return 0;
}