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

- Live demo: https://bicexact.onrender.com/
- Frontend: `frontend/`
- Backend: `src/main.go`
- Estimator executable: `src/estimate_hinge_numbers`
- API endpoint: `POST /api/estimate`

### Prerequisites
- Go 1.26.10+ (use a supported release with current security patches)
- Node.js and npm
- A C++14 compiler with OpenMP support
  - Linux/Docker: `g++`
  - macOS: Homebrew `llvm`, `libomp`, and `eigen`

### Local quick start
Rebuild the estimator executable after updating the source, including security updates. The bundled executable is for macOS ARM64; Linux deployments are built from source by Docker.

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
npm ci
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
npm ci
npm run dev
```

Open the Vite URL printed in the terminal. API requests to `/api` are proxied to `http://localhost:8080`.

### Docker
Build and run the containerized web application.

```bash
docker compose up -d --build
```

Open http://localhost:8080/ in your browser.

The production container runs as UID 10001 and contains only the server, estimator, frontend, and runtime libraries. Compose binds port 8080 to localhost and limits CPU, memory, processes, and temporary storage. Application files are read-only. Use the separate research container below for compilation and Python experiments.

### Request limits and security checks

Both upload endpoints accept a maximum of 64 MiB per request, two concurrent requests, and a 30-second calculation timeout. Uploads must be nonempty `.pdb`, `.ent`, or `.cif` files; chain IDs contain 1–8 ASCII letters or digits. Structures must have matching CA coordinates and at most 10,000 CA residues in each selected chain. Cross-origin browser submissions are rejected. Invalid inputs return 400, oversized requests 413, busy capacity 503, and calculation timeouts 504.

Run the regression tests after building the estimator:

```bash
cd src
go test -race ./...
cd ..
BICEXACT_ESTIMATOR=src/estimate_hinge_numbers python3 -m unittest discover -s tests -v
```

Audit the locked dependencies:

```bash
cd frontend
npm audit
cd ../src
go run golang.org/x/vuln/cmd/govulncheck@latest ./...
cd ..
uv run --with pip-audit scripts/audit_python_dependencies.py
```

The 2026-10-22 audit found no npm or Python vulnerabilities and no vulnerabilities in imported Go packages. Go still reports [GO-2026-5932](https://pkg.go.dev/vuln/GO-2026-5932) at module level: the unmaintained `golang.org/x/crypto/openpgp` package has no fixed version. This application does not import that package; `govulncheck` reports no reachable affected symbols. These checks do not constitute an OS package vulnerability scan.

## Experiments
- Stop the web container before starting the research container, which also uses the name `bel`. The research environment includes compilers, Python 3.12, and Poetry 2.2.1 and mounts `src/` for experiment outputs.
```bash
docker compose down
docker compose -f docker-compose.experiments.yml up -d --build
```
- If you want to work inside the docker container, execute the following command
```bash
docker exec -it bel bash
```
### Install libraries for Python
- inside the docker container
```bash
poetry install --no-root
```
  - outside the docker container
```bash
docker exec bel poetry install --no-root
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
