#pragma once

#include <cstdint>
#include <vector>
#include <unordered_map>

struct Slot {
    bool occupied;
    uint64_t key;
};

using HistogramFixed = std::vector<size_t>;
using HistogramDynamic = std::unordered_map<size_t, size_t>;