#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

uint64_t hash_uint64(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

struct ProbeHistogram {
    // bins[i] counts probe length == i; bins.back() is an overflow bucket
    std::vector<uint64_t> bins;
    explicit ProbeHistogram(std::size_t max_bin = 512) : bins(max_bin + 1, 0) {}
    void clear() { std::fill(bins.begin(), bins.end(), 0); }
    void record(uint64_t probes) {
        if (probes >= bins.size() - 1) {
            bins.back()++;
        } else {
            bins[probes]++;
        }
    }
    uint64_t total() const {
        uint64_t t = 0;
        for (auto c : bins) t += c;
        return t;
    }
    uint64_t percentile(double p) const {  // p in [0,1]
        uint64_t t = total();
        if (t == 0) return 0;
        uint64_t target = static_cast<uint64_t>(p * t);
        uint64_t run = 0;
        for (std::size_t i = 0; i < bins.size(); ++i) {
            run += bins[i];
            if (run >= target) return static_cast<uint64_t>(i);
        }
        return static_cast<uint64_t>(bins.size() - 1);
    }
    const std::vector<uint64_t>& data() const { return bins; }
};

struct Bucket {
    bool occupied = false;
    uint64_t key = 0;
    int32_t probe_distance = 0;  // kept for identical layout across tables
};

class LinearProbingHashTable {
public:
    explicit LinearProbingHashTable(std::size_t capacity)
        : buckets_(capacity, Bucket{}), mask_(capacity - 1) {}

    void reset_all_metrics() { reset_insert_metrics(); reset_lookup_metrics(); }
    void reset_insert_metrics() {
        total_insert_probes_ = 0;
        total_insert_ops_ = 0;
        max_insert_probes_ = 0;
        insert_hist_.clear();
    }
    void reset_lookup_metrics() {
        total_lookup_probes_ = 0;
        total_lookup_ops_ = 0;
        max_lookup_probes_ = 0;
        lookup_hist_.clear();
    }

    void insert(uint64_t key) {
        std::size_t index = hash_uint64(key) & mask_;
        uint64_t probes = 0;

        while (true) {
            Bucket& bucket = buckets_[index];
            ++probes;

            if (!bucket.occupied) {
                bucket.occupied = true;
                bucket.key = key;
                bucket.probe_distance = 0;
                break;
            }
            if (bucket.key == key) break;  // already present
            index = (index + 1) & mask_;
        }

        total_insert_probes_ += probes;
        ++total_insert_ops_;
        if (probes > max_insert_probes_) max_insert_probes_ = probes;
        insert_hist_.record(probes);
    }

    bool lookup(uint64_t key) {
        std::size_t index = hash_uint64(key) & mask_;
        uint64_t probes = 0;
        bool found = false;

        while (true) {
            const Bucket& bucket = buckets_[index];
            ++probes;

            if (!bucket.occupied) break;
            if (bucket.key == key) { found = true; break; }
            index = (index + 1) & mask_;
        }

        total_lookup_probes_ += probes;
        ++total_lookup_ops_;
        if (probes > max_lookup_probes_) max_lookup_probes_ = probes;
        lookup_hist_.record(probes);
        return found;
    }

    double avg_insert_probes() const {
        return total_insert_ops_ == 0 ? 0.0
               : static_cast<double>(total_insert_probes_) /
                     static_cast<double>(total_insert_ops_);
    }
    double avg_lookup_probes() const {
        return total_lookup_ops_ == 0 ? 0.0
               : static_cast<double>(total_lookup_probes_) /
                     static_cast<double>(total_lookup_ops_);
    }

    uint64_t max_insert_probes() const { return max_insert_probes_; }
    uint64_t max_lookup_probes() const { return max_lookup_probes_; }
    const ProbeHistogram& insert_hist() const { return insert_hist_; }
    const ProbeHistogram& lookup_hist() const { return lookup_hist_; }

private:
    std::vector<Bucket> buckets_;
    std::size_t mask_;

    uint64_t total_insert_probes_ = 0;
    uint64_t total_insert_ops_ = 0;
    uint64_t max_insert_probes_ = 0;
    ProbeHistogram insert_hist_{};

    uint64_t total_lookup_probes_ = 0;
    uint64_t total_lookup_ops_ = 0;
    uint64_t max_lookup_probes_ = 0;
    ProbeHistogram lookup_hist_{};
};

class RobinHoodHashTable {
public:
    explicit RobinHoodHashTable(std::size_t capacity)
        : buckets_(capacity, Bucket{}), mask_(capacity - 1) {}

    void reset_all_metrics() { reset_insert_metrics(); reset_lookup_metrics(); }
    void reset_insert_metrics() {
        total_insert_probes_ = 0;
        total_insert_ops_ = 0;
        max_insert_probes_ = 0;
        insert_hist_.clear();
    }
    void reset_lookup_metrics() {
        total_lookup_probes_ = 0;
        total_lookup_ops_ = 0;
        max_lookup_probes_ = 0;
        lookup_hist_.clear();
    }

    void insert(uint64_t key) {
        std::size_t index = hash_uint64(key) & mask_;
        int32_t probe_distance = 0;
        uint64_t probes = 0;

        while (true) {
            Bucket& bucket = buckets_[index];
            ++probes;

            if (!bucket.occupied) {
                bucket.occupied = true;
                bucket.key = key;
                bucket.probe_distance = probe_distance;
                break;
            }

            if (bucket.key == key) {
                break;  // already present
            }

            if (bucket.probe_distance < probe_distance) {
                std::swap(key, bucket.key);
                std::swap(probe_distance, bucket.probe_distance);
            }

            index = (index + 1) & mask_;
            ++probe_distance;
        }

        total_insert_probes_ += probes;
        ++total_insert_ops_;
        if (probes > max_insert_probes_) max_insert_probes_ = probes;
        insert_hist_.record(probes);
    }

    bool lookup(uint64_t key) {
        std::size_t index = hash_uint64(key) & mask_;
        int32_t probe_distance = 0;
        uint64_t probes = 0;
        bool found = false;

        while (true) {
            const Bucket& bucket = buckets_[index];
            ++probes;

            if (!bucket.occupied) break;
            if (bucket.probe_distance < probe_distance) break;  // would be earlier
            if (bucket.key == key) { found = true; break; }

            index = (index + 1) & mask_;
            ++probe_distance;
        }

        total_lookup_probes_ += probes;
        ++total_lookup_ops_;
        if (probes > max_lookup_probes_) max_lookup_probes_ = probes;
        lookup_hist_.record(probes);
        return found;
    }

    double avg_insert_probes() const {
        return total_insert_ops_ == 0 ? 0.0
               : static_cast<double>(total_insert_probes_) /
                     static_cast<double>(total_insert_ops_);
    }
    double avg_lookup_probes() const {
        return total_lookup_ops_ == 0 ? 0.0
               : static_cast<double>(total_lookup_probes_) /
                     static_cast<double>(total_lookup_ops_);
    }

    uint64_t max_insert_probes() const { return max_insert_probes_; }
    uint64_t max_lookup_probes() const { return max_lookup_probes_; }
    const ProbeHistogram& insert_hist() const { return insert_hist_; }
    const ProbeHistogram& lookup_hist() const { return lookup_hist_; }

private:
    std::vector<Bucket> buckets_;
    std::size_t mask_;

    uint64_t total_insert_probes_ = 0;
    uint64_t total_insert_ops_ = 0;
    uint64_t max_insert_probes_ = 0;
    ProbeHistogram insert_hist_{};

    uint64_t total_lookup_probes_ = 0;
    uint64_t total_lookup_ops_ = 0;
    uint64_t max_lookup_probes_ = 0;
    ProbeHistogram lookup_hist_{};
};

std::vector<uint64_t> generate_unique_keys(std::size_t count, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<uint64_t> dist;

    std::unordered_set<uint64_t> seen;
    seen.reserve(count * 2);

    std::vector<uint64_t> keys;
    keys.reserve(count);

    while (keys.size() < count) {
        uint64_t candidate = dist(rng);
        if (seen.insert(candidate).second) keys.push_back(candidate);
    }
    return keys;
}

std::vector<uint64_t> generate_absent_keys(std::size_t count,
                                           uint64_t seed,
                                           const std::unordered_set<uint64_t>& forbidden) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<uint64_t> dist;

    std::unordered_set<uint64_t> seen;
    seen.reserve(count * 2);

    std::vector<uint64_t> keys;
    keys.reserve(count);

    while (keys.size() < count) {
        uint64_t candidate = dist(rng);
        if (forbidden.find(candidate) == forbidden.end() &&
            seen.insert(candidate).second) {
            keys.push_back(candidate);
        }
    }
    return keys;
}

struct DualPrinter {
    std::ostream& term;
    std::ostream& file;
    void row(int capacity_log2,
             const std::string& table_type,
             double load_factor,
             int repetition,
             const std::string& phase,
             std::size_t num_ops,
             double time_seconds,
             double avg_probes,
             uint64_t max_probes) {
        auto emit = [&](std::ostream& os) {
            os << capacity_log2 << ','
               << table_type << ','
               << load_factor << ','
               << repetition << ','
               << phase << ','
               << num_ops << ','
               << time_seconds << ','
               << avg_probes << ','
               << max_probes << '\n';
        };
        emit(term);
        emit(file);
    }
};

void write_histogram(std::ostream& os,
                     int capacity_log2,
                     const std::string& table_type,
                     double load_factor,
                     int repetition,
                     const std::string& phase,
                     const ProbeHistogram& hist) {
    uint64_t p50 = hist.percentile(0.50);
    uint64_t p90 = hist.percentile(0.90);
    uint64_t p99 = hist.percentile(0.99);
    const auto& bins = hist.data();
    for (std::size_t i = 0; i < bins.size(); ++i) {
        uint64_t count = bins[i];
        if (count == 0) continue;
        os << capacity_log2 << ','
           << table_type << ','
           << load_factor << ','
           << repetition << ','
           << phase << ','
           << i << ','
           << count << ','
           << p50 << ','
           << p90 << ','
           << p99 << '\n';
    }
}

}  // namespace

int main() {
    const std::vector<int> capacity_log2_values = {18, 20, 22};
    const std::vector<double> load_factors = {0.50, 0.70, 0.80, 0.85, 0.90, 0.93, 0.95, 0.97};
    const int NUM_REPETITIONS = 30;
    const uint64_t base_seed = 0x12345678abcdef00ULL;

    std::ofstream csv("results.csv", std::ios::trunc);
    if (!csv) {
        std::cerr << "Failed to open results.csv for writing\n";
        return 1;
    }
    std::ofstream histcsv("probe_hist.csv", std::ios::trunc);
    if (!histcsv) {
        std::cerr << "Failed to open probe_hist.csv for writing\n";
        return 1;
    }

    std::cout << std::fixed << std::setprecision(6);
    csv << std::fixed << std::setprecision(6);
    histcsv << std::fixed << std::setprecision(6);

    // CSV headers
    std::cout << "capacity_log2,table_type,load_factor,repetition,phase,num_ops,time_seconds,avg_probes,max_probes\n";
    csv      << "capacity_log2,table_type,load_factor,repetition,phase,num_ops,time_seconds,avg_probes,max_probes\n";
    histcsv  << "capacity_log2,table_type,load_factor,repetition,phase,probe_length,count,p50,p90,p99\n";

    DualPrinter out{std::cout, csv};

    for (int capacity_log2 : capacity_log2_values) {
        std::size_t capacity = 1ull << capacity_log2;

        for (std::size_t lf_index = 0; lf_index < load_factors.size(); ++lf_index) {
            double load_factor = load_factors[lf_index];
            std::size_t num_keys = static_cast<std::size_t>(load_factor * static_cast<double>(capacity));

            for (int repetition = 1; repetition <= NUM_REPETITIONS; ++repetition) {
                uint64_t insert_seed = base_seed
                    + static_cast<uint64_t>(capacity_log2) * 100000ULL
                    + static_cast<uint64_t>(lf_index) * 1000ULL
                    + static_cast<uint64_t>(repetition);

                std::vector<uint64_t> insert_keys = generate_unique_keys(num_keys, insert_seed);
                std::unordered_set<uint64_t> inserted_set(insert_keys.begin(), insert_keys.end());

                uint64_t miss_seed = insert_seed + 0x9e3779b97f4a7c15ULL;
                std::vector<uint64_t> missing_keys =
                    generate_absent_keys(num_keys, miss_seed, inserted_set);

                // Linear probing
                {
                    LinearProbingHashTable table(capacity);

                    // Insert phase
                    table.reset_all_metrics();
                    auto start = std::chrono::high_resolution_clock::now();
                    for (uint64_t key : insert_keys) table.insert(key);
                    auto end = std::chrono::high_resolution_clock::now();
                    double time_seconds = std::chrono::duration<double>(end - start).count();
                    out.row(capacity_log2, "linear", load_factor, repetition, "insert",
                            num_keys, time_seconds, table.avg_insert_probes(), table.max_insert_probes());
                    write_histogram(histcsv, capacity_log2, "linear", load_factor, repetition,
                                    "insert", table.insert_hist());

                    // Lookup success phase
                    table.reset_lookup_metrics();
                    start = std::chrono::high_resolution_clock::now();
                    for (uint64_t key : insert_keys) table.lookup(key);
                    end = std::chrono::high_resolution_clock::now();
                    time_seconds = std::chrono::duration<double>(end - start).count();
                    out.row(capacity_log2, "linear", load_factor, repetition, "lookup_success",
                            num_keys, time_seconds, table.avg_lookup_probes(), table.max_lookup_probes());
                    write_histogram(histcsv, capacity_log2, "linear", load_factor, repetition,
                                    "lookup_success", table.lookup_hist());

                    // Lookup fail phase
                    table.reset_lookup_metrics();
                    start = std::chrono::high_resolution_clock::now();
                    for (uint64_t key : missing_keys) table.lookup(key);
                    end = std::chrono::high_resolution_clock::now();
                    time_seconds = std::chrono::duration<double>(end - start).count();
                    out.row(capacity_log2, "linear", load_factor, repetition, "lookup_fail",
                            num_keys, time_seconds, table.avg_lookup_probes(), table.max_lookup_probes());
                    write_histogram(histcsv, capacity_log2, "linear", load_factor, repetition,
                                    "lookup_fail", table.lookup_hist());
                }

                // Robin Hood hashing
                {
                    RobinHoodHashTable table(capacity);

                    // Insert phase
                    table.reset_all_metrics();
                    auto start = std::chrono::high_resolution_clock::now();
                    for (uint64_t key : insert_keys) table.insert(key);
                    auto end = std::chrono::high_resolution_clock::now();
                    double time_seconds = std::chrono::duration<double>(end - start).count();
                    out.row(capacity_log2, "robinhood", load_factor, repetition, "insert",
                            num_keys, time_seconds, table.avg_insert_probes(), table.max_insert_probes());
                    write_histogram(histcsv, capacity_log2, "robinhood", load_factor, repetition,
                                    "insert", table.insert_hist());

                    // Lookup success phase
                    table.reset_lookup_metrics();
                    start = std::chrono::high_resolution_clock::now();
                    for (uint64_t key : insert_keys) table.lookup(key);
                    end = std::chrono::high_resolution_clock::now();
                    time_seconds = std::chrono::duration<double>(end - start).count();
                    out.row(capacity_log2, "robinhood", load_factor, repetition, "lookup_success",
                            num_keys, time_seconds, table.avg_lookup_probes(), table.max_lookup_probes());
                    write_histogram(histcsv, capacity_log2, "robinhood", load_factor, repetition,
                                    "lookup_success", table.lookup_hist());

                    // Lookup fail phase
                    table.reset_lookup_metrics();
                    start = std::chrono::high_resolution_clock::now();
                    for (uint64_t key : missing_keys) table.lookup(key);
                    end = std::chrono::high_resolution_clock::now();
                    time_seconds = std::chrono::duration<double>(end - start).count();
                    out.row(capacity_log2, "robinhood", load_factor, repetition, "lookup_fail",
                            num_keys, time_seconds, table.avg_lookup_probes(), table.max_lookup_probes());
                    write_histogram(histcsv, capacity_log2, "robinhood", load_factor, repetition,
                                    "lookup_fail", table.lookup_hist());
                }
            }
        }
    }

    csv.flush();
    histcsv.flush();
    return 0;
}
