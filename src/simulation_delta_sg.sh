#!/bin/bash
set -eu

"${CXX:-g++}" simulation_delta_g_parallel_latest.cpp -o simulation_delta_g_parallel_latest -std=c++14 -Wall -Wextra -O3 -Werror -fopenmp
./simulation_delta_g_parallel_latest "output_simulation_file.csv"
