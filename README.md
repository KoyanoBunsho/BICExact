# BIC Exact and BIC LH
- Fast and accurate algorithms for estimating the number of hinges in a protein based on information criteria.
## Input and Output
- Input: two conformations of the same protein
- Output: the number of hinges and the corresponding hinge positions
## Methods
- For details on the methods, please refer to our upcoming paper.
### Notations
- $n$: the number of observed data points
- $H^{(k)} = [h_0, h_1, \dots, h_{k+1}]$: a list of $k + 2$ distinct positive integers
- $R^{(k)} = [\boldsymbol{R}_1, \dots, \boldsymbol{R}_{k+1}]$: a list of $k+1$ rotation matrices
- $V^{(k)} = [\boldsymbol{v}_1, \dots, \boldsymbol{v}_{k+1}]$: a list of $k+1$ translation vectors
- $\boldsymbol{\epsilon}_i = \begin{pmatrix}
    \epsilon_{i, x}\\ \epsilon_{i, y}\\ \epsilon_{i, z}
\end{pmatrix} \in \mathbb{R}^3$ be the $i$-th error term
### $k$-hinge model

$$
\boldsymbol{y}_i = f_{R^{(k)}, V^{(k)}, H^{(k)}}(\boldsymbol{x}_i, i) + \boldsymbol{\epsilon}_i
$$

### Hinge number estimation

$$
k^{(A, B, c)} = argmin_{0 \leq k \leq n-1}\min_{R^{(k)}, V^{(k)}, H^{(k)}}\sum_{r=1}^{k+1}\{\sum_{i=h_{r - 1}}^{h_{r}-1}||\boldsymbol{a}_i - \boldsymbol{R}_r(\boldsymbol{b}_i + \boldsymbol{v}_r)||^2 + c\}
$$

When $c = \sigma^2q\log n$, $k^{(A, B, c)}$-hinge model minimizes BIC.
When $c = 2\sigma^2q$, $k^{(A, B, c)}$-hinge model minimizes AIC.

Using the dynamic programming (DP), we can obtain the exact solution in $O(n^2)$-time.

By applying LARSCH algorithm, we can obtain the heuristic solution in $O(n)$-time.

## Experiments
- build the docker container
```bash
docker-compose up -d
```
- If you want to work inside the docker container, execute the following command
```bash
docker exec -it bel bash
```
### Install libraries for Python
- inside the docker container
```bash
poetry install
```
  - outside the docker container
```bash
docker exec bel poetry install
```
### Constructing a simulation dataset
- Please make sure you downloaded all pdb files in the all_pdb directory before you execute the following command
- inside the docker container
```bash
poetry run python make_simulation_dataset.py
```
- outside the docker container
```bash
docker exec bel poetry run python make_simulation_dataset.py
```
### Experiments on the simulation dataset
- inside the docker container
```bash
bash simulation_bic.sh
```
- outside the docker container
```bash
docker exec bel bash simulation_bic.sh
```
### Experiments on the Shibuya 2008 dataset
- inside the docker container
```bash
bash speed_comparison_rmsdh.sh
```
- outside the docker container
```bash
docker exec bel bash speed_comparison_rmsdh.sh
```
### Evaluating the results on the simulation dataset
- inside the docker container
```bash
poetry run python plot_accuracy_f_measure.py
```
- outside the docker container
```bash
docker exec bel poetry run python plot_accuracy_f_measure.py
```
### Evaluating the results on the Shibuya 2008 dataset
- inside the docker container
```bash
poetry run python eval_result_shibuya.py
```
- outside the docker container
```bash
docker exec bel poetry run python eval_result_shibuya.py
```

## Example usage
### CUI
#### When you execute the command outside the docker container
```bash
docker exec bel g++ estimate_hinge_numbers.cpp -o estimate_hinge_numbers -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
docker exec bel ./estimate_hinge_numbers pdb3hvp.pdb pdb4hvp.pdb A A bic exact
docker exec bel ./estimate_hinge_numbers pdb3hvp.pdb pdb4hvp.pdb A A bic lh
```
#### When you execute the command inside the docker container
```bash
g++ estimate_hinge_numbers.cpp -o estimate_hinge_numbers -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
./estimate_hinge_numbers pdb3hvp.pdb pdb4hvp.pdb A A bic exact
./estimate_hinge_numbers pdb3hvp.pdb pdb4hvp.pdb A A bic lh
```
### GUI
- After you build the estimate_hinge_numbers.cpp, start the web server on your local computer
```bash
cd src
go mod download
go run main.go
```
- access the following link: http://localhost:8080/

![Demonstration of the application](demo.gif)
