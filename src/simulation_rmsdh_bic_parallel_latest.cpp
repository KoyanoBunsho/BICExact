//go:build ignore

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("unroll-loops")
#endif

#include "rmsd_io.h"
#include "rmsd_struct.h"
#include "rmsdh_new.h"
#include <sstream>

std::string extractHingeIndices(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "ファイルが開けません: " << filename << std::endl;
    return "";
  }
  std::string line;
  std::getline(file, line);
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string cell;
    std::vector<std::string> row;
    while (std::getline(ss, cell, ',')) {
      row.push_back(cell);
    }
    if (row.size() > 5) {
      file.close();
      return row[4];
    }
  }
  file.close();
  return "";
}

int main(int argc, char **argv) try {
  if (argc != 4) throw std::invalid_argument("Expected sigma, method, and simulation CSV");
  const double sigma = read_sigma(argv[1]);
  bool is_heuristic = false;
  std::string save_method_name;
  std::ofstream myfile;
  int skipped = 0;
  if (std::string(argv[2]) == "aic") {
    save_method_name = "rmsdh_aic";
  } else if (std::string(argv[2]) == "bic") {
    save_method_name = "rmsdh_bic";
  } else if (std::string(argv[2]) == "lh_aic") {
    save_method_name = "fast_rmsdh_aic";
    is_heuristic = true;
  } else if (std::string(argv[2]) == "lh_bic") {
    save_method_name = "fast_rmsdh_bic";
    is_heuristic = true;
  } else {
    throw std::invalid_argument("Unknown estimation method");
  }
  std::string input_file_path = std::string(argv[3]);
  const auto file_triples = read_simulation_files(input_file_path);
  std::string save_name = "rmsdh_result/simulation_parallel_" +
                          save_method_name + "_" + std::to_string(sigma) +
                          input_file_path.substr(input_file_path.find_last_of("/")+1);
  myfile.open(save_name);
  if (!myfile) throw std::runtime_error("Cannot open results file");
  myfile << "p_pdb_id,Residue "
            "length,RMSD,RMSDh,c,actual_hinge_cnt,hinge_cnt,actual_hinge_"
            "indices,hinge_index,exec_time (s)"
         << std::endl;
  const bool use_bic = std::string(argv[2]) == "bic" || std::string(argv[2]) == "lh_bic";
#pragma omp parallel for reduction(+ : skipped)
  for (const auto &triple : file_triples) {
    std::string p_pdb_id =
        std::get<0>(triple).substr(std::get<0>(triple).find_last_of("/") + 1);
    Eigen::MatrixXd p = openMatrixData(std::get<0>(triple));
    Eigen::MatrixXd q = openMatrixData(std::get<1>(triple));
    std::string hinge_file = std::get<2>(triple);
    std::string hingeIndices = extractHingeIndices(hinge_file);
    int total_residue_length = p.cols();
    if (p.cols() == 0 || p.cols() != q.cols()) {
      std::cout << "p length: " << p.cols() << " q length: " << q.cols()
                << std::endl;
      std::cerr << "Skipping empty, invalid, or mismatched coordinate matrices" << std::endl;
      ++skipped;
      continue;
    }
    std::cout << total_residue_length << std::endl;
    std::vector<double> default_weights(total_residue_length, 1.0);
    ConformationPair PQ_pair = MoveToOrigin(p, q, default_weights);
    double rmsd_result = CalcRMSD(PQ_pair.P, PQ_pair.Q, default_weights);
    std::cout << p_pdb_id << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    const double c = 7 * sigma * sigma * (use_bic ? std::log(total_residue_length) : 2.0);
    ProteinRMSDhinge rmsdh_calculator(PQ_pair.P, PQ_pair.Q, c);
    RMSDhHingeCnt rmsdh_hinge_cnt_result;
    if (is_heuristic) {
      rmsdh_hinge_cnt_result = rmsdh_calculator.CalcLinearRMSDh();
    } else {
      rmsdh_hinge_cnt_result = rmsdh_calculator.CalcRMSDh();
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> exec_time_ms = end - start;
    double exec_time_s = exec_time_ms.count() / 1000.0;
    std::cout << exec_time_s << " s" << std::endl;
    double rmsdh_result = rmsdh_hinge_cnt_result.rmsdh_result;
    int hinge_cnt = rmsdh_hinge_cnt_result.hinge_cnt;
    std::vector<int> hinge_index_vec = rmsdh_hinge_cnt_result.hinge_index_vec;
    std::string hinge_index = "";
    for (int i = hinge_index_vec.size() - 1; i >= 0; i--) {
      if (i > 0) {
        hinge_index += (std::to_string(hinge_index_vec[i]) + " : ");
      } else {
        hinge_index += (std::to_string(hinge_index_vec[i]));
      }
    }
#pragma omp critical
    {
      myfile << p_pdb_id << ",";
      myfile << total_residue_length << ",";
      myfile << rmsd_result << ",";
      myfile << rmsdh_result << "," << c << "," << hinge_file << ","
             << hinge_cnt << "," << hingeIndices << "," << hinge_index << ","
             << exec_time_s << std::endl;
    }
  }
  myfile.close();
  return myfile && skipped == 0 ? 0 : 2;
} catch (const std::exception &err) {
  std::cerr << "Invalid experiment input: " << err.what() << '\n';
  return 2;
}
