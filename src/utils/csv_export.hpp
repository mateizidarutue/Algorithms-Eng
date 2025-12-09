#pragma once

#include <string>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <cstdlib> // system()

namespace csvutil {

inline void ensure_output_dir() {
    system("if not exist output mkdir output");
}

inline void export_fixed_histogram(
        const std::string& filename,
        const std::vector<size_t>& hist)
{
    ensure_output_dir();
    std::ofstream out("output/" + filename);

    out << "probes,count\n";
    for (size_t i = 0; i < hist.size(); ++i) {
        out << i << "," << hist[i] << "\n";
    }
}

inline void export_dynamic_histogram(
        const std::string& filename,
        const std::unordered_map<size_t,size_t>& hist)
{
    ensure_output_dir();
    std::ofstream out("output/" + filename);

    out << "probes,count\n";
    for (auto& [dist, count] : hist) {
        out << dist << "," << count << "\n";
    }
}

} // namespace csvutil
