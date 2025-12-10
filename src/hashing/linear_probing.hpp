#pragma once

#include <vector>
#include <iostream>
#include <unordered_map>

#include "hash_function.hpp"
#include "types.hpp"
#include "utils/csv_export.hpp"

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
        // CONFIGURABLE: Histogram depth: Bucket count set to 21, raise to capture more probe lengths
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
        std::cout << "  Avg insert probes: " << average_insert_probes() << "\n";
        std::cout << "  Max insert probes: " << max_probe_insert << "\n";
        std::cout << "  Insert variance: " << variance_insert_probes() << "\n";
        std::cout << "  Avg lookup probes: " << average_lookup_probes() << "\n";
        std::cout << "  Lookup variance: " << variance_lookup_probes() << "\n";
        std::cout << "  Avg lookup HIT probes: " << lookup_hit_avg() << "\n";
        std::cout << "  Avg lookup MISS probes: " << lookup_miss_avg() << "\n";

        std::cout<< "----------------------------" << '\n' << '\n';
    }

    void export_histograms_csv(const std::string& prefix) const {
        csvutil::export_fixed_histogram(prefix + "_insert_fixed.csv",
                                        insert_histogram_fixed);
        csvutil::export_fixed_histogram(prefix + "_lookup_fixed.csv",
                                        lookup_histogram_fixed);

        csvutil::export_dynamic_histogram(prefix + "_insert_dynamic.csv",
                                        insert_histogram_dynamic);
        csvutil::export_dynamic_histogram(prefix + "_lookup_dynamic.csv",
                                        lookup_histogram_dynamic);
        }


private:
    // Represents a slot in the hash table, indicating whether it is occupied and storing the key if it is.
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

    /*
        Fixed histograms are used to capture the distribution of probe lengths
        within a predefined range (e.g., 0-20). This allows efficient tracking
        of probe statistics without dynamically resizing the histogram.
    */
    HistogramFixed insert_histogram_fixed;
    HistogramFixed lookup_histogram_fixed;

    /*
        Dynamic histograms use unordered maps to capture the full range of probe lengths.
        This is useful for analyzing outliers and understanding the complete distribution
        of probe lengths beyond the fixed histogram's limits.
    */
    HistogramDynamic insert_histogram_dynamic;
    HistogramDynamic lookup_histogram_dynamic;
};
