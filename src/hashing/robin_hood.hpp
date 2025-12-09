#pragma once

#include <vector>
#include <iostream>
#include <cmath>

#include "hash_function.hpp"
#include "types.hpp"
#include "utils/csv_export.hpp"

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

        std::cout<< "----------------------------" << '\n' << '\n';
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

    void export_histograms_csv(const std::string& prefix) const {
        csvutil::export_fixed_histogram(prefix + "_insert_fixed.csv",
                                        insert_histogram_fixed);
        csvutil::export_fixed_histogram(prefix + "_lookup_fixed.csv",
                                        lookup_histogram_fixed);
        csvutil::export_fixed_histogram(prefix + "_final_distance_fixed.csv",
                                        final_distance_histogram_fixed);

        csvutil::export_dynamic_histogram(prefix + "_insert_dynamic.csv",
                                        insert_histogram_dynamic);
        csvutil::export_dynamic_histogram(prefix + "_lookup_dynamic.csv",
                                        lookup_histogram_dynamic);
        csvutil::export_dynamic_histogram(prefix + "_final_distance_dynamic.csv",
                                        final_distance_histogram_dynamic);
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

    /*
        Fixed histograms are used to capture the distribution of probe lengths
        within a predefined range (e.g., 0-20). This allows efficient tracking
        of probe statistics without dynamically resizing the histogram.
    */
    HistogramFixed insert_histogram_fixed;
    HistogramFixed lookup_histogram_fixed;

    HistogramFixed final_distance_histogram_fixed = std::vector<size_t>(21, 0);
    HistogramDynamic final_distance_histogram_dynamic;

    /*
        Dynamic histograms use unordered maps to capture the full range of probe lengths.
        This is useful for analyzing outliers and understanding the complete distribution
        of probe lengths beyond the fixed histogram's limits.
    */
    HistogramDynamic insert_histogram_dynamic;
    HistogramDynamic lookup_histogram_dynamic;
};