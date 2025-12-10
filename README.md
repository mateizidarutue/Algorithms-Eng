# 2IAB00-Lab-Project

## Configurable Parameters

All configurable experiment variables are marked with:

```cpp
// CONFIGURABLE:

Use CTRL+SHIFT+F and search for CONFIGURABLE to edit them

How to Run

Compile and run:

g++ -std=c++17 -O3 -I./src src/main.cpp -o benchmark_app
./benchmark_app

The program will:
- Generate random keys
- Run Linear Probing and Robin Hood benchmarks
- Print timing + probe statistics
- Export histogram CSV files in output/

CSV files are generated for both hash tables:
- *_insert_fixed.csv
- *_lookup_fixed.csv
- *_insert_dynamic.csv
- *_lookup_dynamic.csv

These files show how many operations required each probe length. To be plotted to visualize clustering and performance differences.

