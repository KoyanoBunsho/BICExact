# BIC Exact and BIC LH
- Fast and accurate algorithms for estimating the number of hinges in a protein based on information criteria.
- $O(n^2)$-time exact algorithm, where $n$ is the protein length
- $O(n)$-time heuristic algorithm, where $n$ is the protein length

## Input and Output
- Input: two conformations of the same protein
- Output: the number of hinges and the corresponding hinge positions

## Methods
- For details on the methods, please refer to our upcoming paper.

## Web application
The repository includes a React + TypeScript frontend and a Go backend.

- Frontend: `frontend/`
- Backend: `src/main.go`
- Estimator executable: `src/estimate_hinge_numbers`
- API endpoint: `POST /api/estimate`

### Prerequisites
- Go 1.20+
- Node.js and npm
- A C++14 compiler with OpenMP support
  - Linux/Docker: `g++`
  - macOS: Homebrew `llvm`, `libomp`, and `eigen`

### Local quick start
Build the estimator executable first.

On Linux or inside the Docker container:

```bash
cd src
g++ estimate_hinge_numbers.cpp -o estimate_hinge_numbers -std=c++14 -lstdc++fs -Wall -Wextra -O3 -Werror -fopenmp
```

On macOS, Apple clang does not support `-fopenmp`. Install Homebrew LLVM, OpenMP, and Eigen, then compile with Homebrew `clang++`.

```bash
brew install llvm libomp eigen
cd src
$(brew --prefix llvm)/bin/clang++ estimate_hinge_numbers.cpp -o estimate_hinge_numbers -std=c++14 -Wall -Wextra -O3 -Werror -fopenmp -I$(brew --prefix libomp)/include -I$(brew --prefix eigen)/include -L$(brew --prefix libomp)/lib -lomp
```

Install and build the React frontend.

```bash
cd ../frontend
npm install
npm run build
```

Start the Go server.

```bash
cd ../src
go mod download
go run main.go
```

Open http://localhost:8080/ in your browser.

### Frontend development
For UI development, keep the Go server running on port `8080` and start Vite in another terminal.

```bash
cd frontend
npm install
npm run dev
```

Open the Vite URL printed in the terminal. API requests to `/api` are proxied to `http://localhost:8080`.

### Docker
Build and run the containerized web application.

```bash
docker-compose up -d
```

Open http://localhost:8080/ in your browser.

To work inside the container:

```bash
docker exec -it bel bash
```

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
### Hinge estimation examples
```bash
bash hinge_estimation_example.sh
```

## Example usage
### CUI
#### When you execute the command outside the docker container
```bash
docker exec bel g++ estimate_hinge_numbers.cpp -o estimate_hinge_numbers -std=c++14 -lstdc++fs -Wall -Wextra -O3 -Werror -fopenmp
docker exec bel ./estimate_hinge_numbers pdb3hvp.pdb pdb4hvp.pdb A A bic exact
docker exec bel ./estimate_hinge_numbers pdb3hvp.pdb pdb4hvp.pdb A A bic lh
```
#### When you execute the command inside the docker container
```bash
g++ estimate_hinge_numbers.cpp -o estimate_hinge_numbers -std=c++14 -lstdc++fs -Wall -Wextra -O3 -Werror -fopenmp
./estimate_hinge_numbers pdb3hvp.pdb pdb4hvp.pdb A A bic exact
./estimate_hinge_numbers pdb3hvp.pdb pdb4hvp.pdb A A bic lh
```
#### When you execute the command on macOS
```bash
brew install llvm libomp eigen
$(brew --prefix llvm)/bin/clang++ estimate_hinge_numbers.cpp -o estimate_hinge_numbers -std=c++14 -Wall -Wextra -O3 -Werror -fopenmp -I$(brew --prefix libomp)/include -I$(brew --prefix eigen)/include -L$(brew --prefix libomp)/lib -lomp
./estimate_hinge_numbers pdb3hvp.pdb pdb4hvp.pdb A A bic exact
./estimate_hinge_numbers pdb3hvp.pdb pdb4hvp.pdb A A bic lh
```
### GUI
See [Web application](#web-application).

![Demonstration of the application](demo.gif)
