//go:build ignore

#include "rmsdh_new.h"

int main(int argc, char **argv) try {
  if (argc != 2)
    throw std::invalid_argument("Expected a simulation CSV");

  const std::string input_path = argv[1];
  const auto files = read_simulation_files(input_path);
  const auto filename = input_path.substr(input_path.find_last_of('/') + 1);
  std::ofstream results("rmsdh_result/simulation_delta_g_" + filename);
  if (!results)
    throw std::runtime_error("Cannot open results file");
  results << "p_pdb_id,Residue length,delta_g,exec_time (s)\n";

  int skipped = 0;
  for (const auto &file : files) {
    const auto &p_path = std::get<0>(file);
    const auto &q_path = std::get<1>(file);
    const Eigen::MatrixXd p = openMatrixData(p_path);
    const Eigen::MatrixXd q = openMatrixData(q_path);
    if (p.cols() == 0 || p.cols() != q.cols()) {
      std::cerr << "Skipping invalid or mismatched coordinate matrices: "
                << p_path << ", " << q_path << '\n';
      ++skipped;
      continue;
    }

    const std::vector<double> weights(p.cols(), 1.0);
    const auto pair = MoveToOrigin(p, q, weights);
    const auto pdb_id = p_path.substr(p_path.find_last_of('/') + 1);
    std::cout << p.cols() << '\n' << pdb_id << std::endl;
    const auto start = std::chrono::steady_clock::now();
    ProteinRMSDhinge calculator(pair.P, pair.Q, 100);
    const double delta_g = calculator.calcDeltaG();
    const std::chrono::duration<double> elapsed =
        std::chrono::steady_clock::now() - start;
    std::cout << elapsed.count() << " s\n";
    results << pdb_id << ',' << p.cols() << ',' << delta_g << ','
            << elapsed.count() << '\n';
  }
  results.close();
  return results && skipped == 0 ? 0 : 2;
} catch (const std::exception &err) {
  std::cerr << "Invalid experiment input: " << err.what() << '\n';
  return 2;
}
