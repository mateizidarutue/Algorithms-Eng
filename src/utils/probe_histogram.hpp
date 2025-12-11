#pragma once

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <map>

class ProbeHistogram {
private:
    std::map<size_t, uint64_t> probe_counts;  // probe_count -> frequency
    size_t total_operations = 0;
    uint64_t sum_probes = 0;
    uint64_t sum_sq_probes = 0;
    size_t max_probes = 0;

public:
    ProbeHistogram() = default;

    void record_probes(size_t probes) {
        probe_counts[probes]++;
        total_operations++;
        sum_probes += probes;
        sum_sq_probes += probes * probes;
        max_probes = std::max(max_probes, probes);
    }

    double get_average() const {
        if (total_operations == 0) return 0.0;
        return static_cast<double>(sum_probes) / static_cast<double>(total_operations);
    }

    double get_variance() const {
        if (total_operations == 0) return 0.0;
        double avg = get_average();
        double mean_sq = static_cast<double>(sum_sq_probes) / static_cast<double>(total_operations);
        return mean_sq - (avg * avg);
    }

    double get_std_dev() const {
        return std::sqrt(get_variance());
    }

    size_t get_max() const {
        return max_probes;
    }

    size_t get_total_operations() const {
        return total_operations;
    }

    double get_percentile(double p) const {
        if (total_operations == 0 || p < 0.0 || p > 1.0) return 0.0;
        
        uint64_t target_count = static_cast<uint64_t>(p * total_operations);
        uint64_t cumulative = 0;
        
        for (const auto& [probes, count] : probe_counts) {
            cumulative += count;
            if (cumulative >= target_count) {
                return static_cast<double>(probes);
            }
        }
        
        return static_cast<double>(max_probes);
    }

    const std::map<size_t, uint64_t>& get_probe_counts() const {
        return probe_counts;
    }

    void to_csv(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) return;
        
        file << "probes,count\n";
        for (const auto& [probes, count] : probe_counts) {
            file << probes << "," << count << "\n";
        }
    }

    void clear() {
        probe_counts.clear();
        total_operations = 0;
        sum_probes = 0;
        sum_sq_probes = 0;
        max_probes = 0;
    }
};

struct ProbeStats {
    double avg = 0.0;
    double variance = 0.0;
    double std_dev = 0.0;
    size_t max_probes = 0;
    double p50 = 0.0;
    double p90 = 0.0;
    double p95 = 0.0;
    double p99 = 0.0;
    double p999 = 0.0;
    size_t total_ops = 0;

    static ProbeStats from_histogram(const ProbeHistogram& hist) {
        ProbeStats stats;
        stats.avg = hist.get_average();
        stats.variance = hist.get_variance();
        stats.std_dev = hist.get_std_dev();
        stats.max_probes = hist.get_max();
        stats.p50 = hist.get_percentile(0.50);
        stats.p90 = hist.get_percentile(0.90);
        stats.p95 = hist.get_percentile(0.95);
        stats.p99 = hist.get_percentile(0.99);
        stats.p999 = hist.get_percentile(0.999);
        stats.total_ops = hist.get_total_operations();
        return stats;
    }
};
