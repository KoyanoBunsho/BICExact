#!/bin/bash

g++ simulation_delta_g_parallel_latest.cpp -o simulation_delta_g_parallel_latest -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp

echo "-----RMSDh-----"
./simulation_delta_g_parallel_latest "output_simulation_file.csv"
