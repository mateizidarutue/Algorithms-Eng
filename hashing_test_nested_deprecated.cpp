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

/* Classic impl of Robin Hood: Insertion:

compute probe distance (current_index - home_index)
while bucket occupied:
    compute occupant distance
    if incoming_distance > occupant_distance:
        swap entries  (steal the spot)
    advance to next bucket
insert item

   Robin Hood: Lookup:

same as LP except you may stop early if:
probe_distance > what key should have

*/
class LinearProbingHashTable
{
public:
    explicit LinearProbingHashTable(size_t size)
        : table(size, {false, 0}), capacity(size),
          total_probes_insert(0), max_probe_insert(0), num_inserts(0),
          total_probes_lookup(0), max_probe_lookup(0), num_lookups(0),
          total_probes_lookup_hit(0), total_probes_lookup_miss(0),
          num_lookup_hit(0), num_lookup_miss(0),
          total_probes_insert_sq(0), total_probes_lookup_sq(0)
    {
        insert_histogram_fixed.assign(21, 0);
        lookup_histogram_fixed.assign(21, 0);
    }

    void insert(uint64_t key) {
        size_t index = custom_hash(key) % capacity;

        size_t probes = 0; // counts how many steps we take for this insertion

        while (table[index].occupied) {
            if (table[index].key == key) {
                return; // Key already exists
            }

            index = (index + 1) % capacity;
            probes++;
        }

        // Update statistics
        total_probes_insert += probes;
        total_probes_insert_sq += (probes * probes);
        if (probes > max_probe_insert) {
            max_probe_insert = probes;
        }

        if (probes < 20) {
            insert_histogram_fixed[probes]++;
        } else {
            insert_histogram_fixed[20]++;
        }

        insert_histogram_dynamic[probes]++;

        num_inserts++;

        table[index] = {true, key};
    }

    // Checks if a given key exists in the hash table.
    bool lookup(uint64_t key) {
        size_t index = custom_hash(key) % capacity;

        size_t probes = 0; // counts the probes for this lookup

        while (table[index].occupied) {
            if (table[index].key == key) {
                total_probes_lookup += probes;
                total_probes_lookup_sq += probes * probes;
                
                if (probes > max_probe_lookup) {
                    max_probe_lookup = probes;
                }

                total_probes_lookup_hit += probes;
                num_lookup_hit++;

                if (probes < 20) {
                    lookup_histogram_fixed[probes]++;
                } else {
                    lookup_histogram_fixed[20]++;
                }

                lookup_histogram_dynamic[probes]++;

                num_lookups++;

                return true;
            }

            index = (index + 1) % capacity;
            probes++;
        }

        // Key not found
        total_probes_lookup += probes;
        total_probes_lookup_sq += probes * probes;

        if (probes > max_probe_lookup) {
            max_probe_lookup = probes;
        }

        total_probes_lookup_miss += probes;
        num_lookup_miss++;

        if (probes < 20) {
            lookup_histogram_fixed[probes]++;
        } else {
            lookup_histogram_fixed[20]++;
        }

        lookup_histogram_dynamic[probes]++;

        num_lookups++;

        return false;
    }

    // Helper functions to compute metrics
    double average_insert_probes() const {
        return num_inserts == 0 ? 0.0 : (double)total_probes_insert / num_inserts;
    }

    double average_lookup_probes() const {
        return num_lookups == 0 ? 0.0 : (double)total_probes_lookup / num_lookups;
    }

    double lookup_hit_avg() const {
        return num_lookup_hit == 0 ? 0.0 : (double)total_probes_lookup_hit / num_lookup_hit;
    }

    double lookup_miss_avg() const {
        return num_lookup_miss == 0 ? 0.0 : (double)total_probes_lookup_miss / num_lookup_miss;
    }

    double variance_insert_probes() const {
        if (num_inserts == 0) return 0.0;
        double mean = average_insert_probes();
        double mean_sq = (double)total_probes_insert_sq / num_inserts;
        return mean_sq - mean * mean;
    }

    double variance_lookup_probes() const {
        if (num_lookups == 0) return 0.0;
        double mean = average_lookup_probes();
        double mean_sq = (double)total_probes_lookup_sq / num_lookups;
        return mean_sq - mean * mean;
    }

    void print_probe_stats() const {
        std::cout << "\n[Linear Probing Stats]\n";
        std::cout << "Avg insert probes: " << average_insert_probes() << "\n";
        std::cout << "Max insert probes: " << max_probe_insert << "\n";
        std::cout << "Insert variance: " << variance_insert_probes() << "\n";
        std::cout << "Avg lookup probes: " << average_lookup_probes() << "\n";
        std::cout << "Lookup variance: " << variance_lookup_probes() << "\n";
        std::cout << "Avg lookup HIT probes: " << lookup_hit_avg() << "\n";
        std::cout << "Avg lookup MISS probes: " << lookup_miss_avg() << "\n";
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

    // Probe statistics
    size_t total_probes_insert;
    size_t max_probe_insert;
    size_t num_inserts;

    size_t total_probes_lookup;
    size_t max_probe_lookup;
    size_t num_lookups;

    size_t total_probes_lookup_hit;
    size_t total_probes_lookup_miss;
    size_t num_lookup_hit;
    size_t num_lookup_miss;

    size_t total_probes_insert_sq;
    size_t total_probes_lookup_sq;

    std::vector<size_t> insert_histogram_fixed;
    std::vector<size_t> lookup_histogram_fixed;

    std::unordered_map<size_t, size_t> insert_histogram_dynamic;
    std::unordered_map<size_t, size_t> lookup_histogram_dynamic;
};

// Classic Robin Hood Hash Table Implementation
class RobinHoodHashTable
{
public:
    explicit RobinHoodHashTable(size_t size)
        : table(size, {false, 0}), capacity(size),
          total_probes_insert(0), max_probe_insert(0), num_inserts(0),
          total_probes_lookup(0), max_probe_lookup(0), num_lookups(0),
          total_probes_lookup_hit(0), total_probes_lookup_miss(0),
          num_lookup_hit(0), num_lookup_miss(0),
          total_probes_insert_sq(0), total_probes_lookup_sq(0)
    {
        insert_histogram_fixed.assign(21, 0);
        lookup_histogram_fixed.assign(21, 0);
    }

    void insert(uint64_t key) {
        size_t index = custom_hash(key) % capacity;

        size_t probe_distance = 0; // how far the current element was displaced
        size_t probes_this_insert = 0; // total probes for this insertion

        Slot current{true, key};  // element that we try to place

        while (true) {
            if (!table[index].occupied) {
                // Found empty spot, place the element here
                table[index] = current;

                size_t final_home = custom_hash(current.key) % capacity;
                size_t final_distance = (index + capacity - final_home) % capacity;

                total_probes_insert += probes_this_insert;
                total_probes_insert_sq += probes_this_insert * probes_this_insert;

                if (probes_this_insert > max_probe_insert) {
                    max_probe_insert = probes_this_insert;
                }

                if (probes_this_insert < 20) {
                    insert_histogram_fixed[probes_this_insert]++;
                } else {
                    insert_histogram_fixed[20]++;
                }
                insert_histogram_dynamic[probes_this_insert]++;

                final_distance_histogram_dynamic[final_distance]++;
                if (final_distance < 20) {
                    final_distance_histogram_fixed[final_distance]++;
                } else {
                    final_distance_histogram_fixed[20]++;
                }

                num_inserts++;

                return;
            }

            // If the key is already present, ignore
            if (table[index].key == key) {
                return;
            }

            // Compute probe distance for the occupant
            size_t occupant_home = custom_hash(table[index].key) % capacity;
            size_t occupant_probe_distance = (index + capacity - occupant_home) % capacity;

            // Robin Hood: steal spot if current element is larger("more deserving")
            if (probe_distance > occupant_probe_distance) {
                std::swap(current, table[index]);
            }

            // Continue probing: move index forward
            index = (index + 1) % capacity;
            probe_distance++;
            probes_this_insert++;
        }
    }

    bool lookup(uint64_t key) {
        size_t index = custom_hash(key) % capacity;

        size_t probe_distance = 0;
        size_t probes_this_lookup = 0;

        while (true) {
            if (!table[index].occupied) {
                total_probes_lookup += probes_this_lookup;
                total_probes_lookup_sq += probes_this_lookup * probes_this_lookup;

                total_probes_lookup_miss += probes_this_lookup;
                num_lookup_miss++;

                if (probes_this_lookup > max_probe_lookup) {
                    max_probe_lookup = probes_this_lookup;
                }

                if (probes_this_lookup < 20) {
                    lookup_histogram_fixed[probes_this_lookup]++;
                } else {
                    lookup_histogram_fixed[20]++;
                }

                lookup_histogram_dynamic[probes_this_lookup]++;

                num_lookups++;
                return false; // Empty spot hit, key not present
            }

            if (table[index].key == key) {
                total_probes_lookup += probes_this_lookup;
                total_probes_lookup_sq += probes_this_lookup * probes_this_lookup;

                total_probes_lookup_hit += probes_this_lookup;
                num_lookup_hit++;

                if (probes_this_lookup > max_probe_lookup) {
                    max_probe_lookup = probes_this_lookup;
                }

                if (probes_this_lookup < 20) {
                    lookup_histogram_fixed[probes_this_lookup]++;
                } else {
                    lookup_histogram_fixed[20]++;
                }

                lookup_histogram_dynamic[probes_this_lookup]++;

                num_lookups++;
                return true;
            }

            // Compute probe distance for the occupant
            size_t occupant_home = custom_hash(table[index].key) % capacity;
            size_t occupant_probe_distance = (index + capacity - occupant_home) % capacity;

            // Early termination: if probe distance exceeds occupant's, key cannot be present
            if (probe_distance > occupant_probe_distance) {
                total_probes_lookup += probes_this_lookup;
                total_probes_lookup_sq += probes_this_lookup * probes_this_lookup;

                total_probes_lookup_miss += probes_this_lookup;
                num_lookup_miss++;

                if (probes_this_lookup > max_probe_lookup) {
                    max_probe_lookup = probes_this_lookup;
                }

                if (probes_this_lookup < 20) {
                    lookup_histogram_fixed[probes_this_lookup]++;
                } else {
                    lookup_histogram_fixed[20]++;
                }

                lookup_histogram_dynamic[probes_this_lookup]++;

                num_lookups++;
                return false;
            }

            // Continue probing
            index = (index + 1) % capacity;
            probe_distance++;
            probes_this_lookup++;
        }
    }

    void print_probe_stats() const {
        std::cout << "\n[Robin Hood Stats]\n";
        std::cout << "Avg insert probes: " << average_insert_probes() << "\n";
        std::cout << "Max insert probes: " << max_probe_insert << "\n";
        std::cout << "Insert variance: " << variance_insert_probes() << "\n";
        std::cout << "Avg lookup probes: " << average_lookup_probes() << "\n";
        std::cout << "Lookup variance: " << variance_lookup_probes() << "\n";
        std::cout << "Avg lookup HIT probes: " << lookup_hit_avg() << "\n";
        std::cout << "Avg lookup MISS probes: " << lookup_miss_avg() << "\n";
    }

    double average_insert_probes() const {
        return num_inserts == 0 ? 0.0 : (double)total_probes_insert / num_inserts;
    }

    double average_lookup_probes() const {
        return num_lookups == 0 ? 0.0 : (double)total_probes_lookup / num_lookups;
    }

    double lookup_hit_avg() const {
        return num_lookup_hit == 0 ? 0.0 : (double)total_probes_lookup_hit / num_lookup_hit;
    }

    double lookup_miss_avg() const {
        return num_lookup_miss == 0 ? 0.0 : (double)total_probes_lookup_miss / num_lookup_miss;
    }

    double variance_insert_probes() const {
        if (num_inserts == 0) return 0.0;

        double mean = average_insert_probes(); // = E[X]
        double mean_sq = (double)total_probes_insert_sq / num_inserts; // = E[X^2]

        return mean_sq - mean * mean; // E[X^2] - (E[X])^2
    }

    double variance_lookup_probes() const {
        if (num_lookups == 0) return 0.0;

        double mean = average_lookup_probes(); // = E[X]
        double mean_sq = (double)total_probes_lookup_sq / num_lookups; // = E[X^2]

        return mean_sq - mean * mean; // E[X^2] - (E[X])^2
    }

private:
    struct Slot {
        bool occupied;
        uint64_t key;
    };

    std::vector<Slot> table;
    size_t capacity;

    size_t total_probes_insert;
    size_t max_probe_insert;
    size_t num_inserts;

    size_t total_probes_lookup;
    size_t max_probe_lookup;
    size_t num_lookups;

    size_t total_probes_lookup_hit;
    size_t total_probes_lookup_miss;
    size_t num_lookup_hit;
    size_t num_lookup_miss;

    size_t total_probes_insert_sq;
    size_t total_probes_lookup_sq;

    std::vector<size_t> insert_histogram_fixed;
    std::vector<size_t> lookup_histogram_fixed;

    std::vector<size_t> final_distance_histogram_fixed = std::vector<size_t>(21, 0);
    std::unordered_map<size_t, size_t> final_distance_histogram_dynamic;

    std::unordered_map<size_t, size_t> insert_histogram_dynamic;
    std::unordered_map<size_t, size_t> lookup_histogram_dynamic;
};

// Encapsulates library implementation of Robin Hood hashing: insertion and lookup timing
void benchmark_robin_hood_library(const std::vector<uint64_t>& keys) {
    // Robin Hood Hashing
    tsl::robin_map<uint64_t, bool, decltype(&custom_hash)> robin_map(0, custom_hash);

    // Measure insertion time for Robin Hood Hashing
    auto start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) {
        robin_map[key] = true;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Robin Hood (library) Insert Time: "
              << std::chrono::duration<double>(end - start).count() << " seconds\n";

    // Measure lookup time for Robin Hood Hashing
    start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) {
        volatile bool found = robin_map.find(key) != robin_map.end();
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Robin Hood (library) Lookup Time: "
              << std::chrono::duration<double>(end - start).count() << " seconds\n";

}

// Encapsulates linear probing insertion and lookup timing
void benchmark_linear_probing(const std::vector<uint64_t>& keys, size_t table_capacity) {
    // Linear Probing Hash Table
    LinearProbingHashTable table(table_capacity);

    // Measure insertion time for Linear Probing
    auto start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) {
        table.insert(key);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "  Linear Probing Insert Time: "
              << std::chrono::duration<double>(end - start).count() << " seconds\n";

    // Measure lookup time for Linear Probing
    start = std::chrono::high_resolution_clock::now();
    for (const auto& key : keys) {
        volatile bool found = table.lookup(key);
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "  Linear Probing Lookup Time: "
              << std::chrono::duration<double>(end - start).count() << " seconds\n";

    table.print_probe_stats();
}

// Benchmarks our implementation of classic robin hood hashing
void benchmark_robinhood_custom(const std::vector<uint64_t>& keys, size_t capacity) {
    RobinHoodHashTable table(capacity);

    // Measure insertion time for classic manual implementation of Robin Hood
    auto start = std::chrono::high_resolution_clock::now();
    for (auto& key : keys) {
        table.insert(key);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Custom Robin Hood Insert Time: "
              << std::chrono::duration<double>(end - start).count() << " seconds\n";

    // Measure lookup time for classic manual implementation of Robin Hood
    start = std::chrono::high_resolution_clock::now();
    for (auto& key : keys) {
        volatile bool found = table.lookup(key);
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Custom Robin Hood Lookup Time: "
              << std::chrono::duration<double>(end - start).count() << " seconds\n";

    table.print_probe_stats();
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

    benchmark_robin_hood_library(keys); // tsl::robin_map
    benchmark_linear_probing(keys, num_elements);
    benchmark_robinhood_custom(keys, num_elements);

    return 0;
}