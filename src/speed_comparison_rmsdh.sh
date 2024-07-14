#!/bin/bash

g++ prmsdh_postprocessing_experiment.cpp -o prmsdh_postprocessing_experiment -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
g++ shibuya_auto_rmsdhk.cpp -o shibuya_auto_rmsdhk -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp

echo "-----RMSDh-----"
locations=("shibuya")

model_types=("aic" "bic")
prefixes=("" "lh_")

for location in "${locations[@]}"; do
    for model_type in "${model_types[@]}"; do
        for prefix in "${prefixes[@]}"; do
            ./prmsdh_postprocessing_experiment 1.55 $location ${prefix}${model_type}
        done
    done
done

./shibuya_auto_rmsdhk
