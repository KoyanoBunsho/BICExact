#!/bin/bash

g++ simulation_rmsdh_bic_parallel_latest.cpp -o simulation_rmsdh_bic_parallel_latest -std=c++14 -lstdc++fs -Wall -Wextra -O3 -Werror -fopenmp
g++ simulation_auto_shibuya_parallel_latest.cpp -o simulation_auto_shibuya_parallel_latest -std=c++14 -lstdc++fs -Wall -Wextra -O3 -Werror -fopenmp

echo "-----RMSDh-----"
model_types=("aic" "bic")
prefixes=("" "lh_")
for sigma in 1.55; do
    for model_type in "${model_types[@]}"; do
        for prefix in "${prefixes[@]}"; do
            ./simulation_rmsdh_bic_parallel_latest $sigma ${prefix}${model_type} "output_simulation_file.csv"
        done
    done
done
./simulation_auto_shibuya_parallel_latest
