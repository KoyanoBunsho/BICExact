#pragma GCC target("avx2")
#pragma GCC optimize("unroll-loops")

#include "rmsd_io.h"
#include "rmsd_struct.h"
#include "rmsdh_new.h"
#include <experimental/filesystem>
#include <regex>
#include <sstream>

namespace fs = std::experimental::filesystem;
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

int main(int argc, char **argv) {
  double sigma = 1.002;
  bool is_heuristic = false;
  std::string save_method_name;
  double c = 0;
  if (argc > 1) {
    sigma = std::stod(argv[1]) / sqrt(3);
  }
  std::string coord_path = "coord_csv_simulation/";
  std::ofstream myfile;
  std::string simulation_data_path = "simulation_data/";
  std::string simulation_data_info_path = "simulation_data_info/";
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
  }
  std::string input_file_path = std::string(argv[3]);
  std::ifstream input_file(input_file_path);
  std::string save_name = "rmsdh_result/simulation_parallel_" +
                          save_method_name + "_" + std::to_string(sigma) +
                          input_file_path;
  myfile.open(save_name);
  myfile << "p_pdb_id,Residue "
            "length,RMSD,RMSDh,c,actual_hinge_cnt,hinge_cnt,actual_hinge_"
            "indices,hinge_index,exec_time (s)"
         << std::endl;
  std::string line;
  std::getline(input_file, line);
  std::vector<std::tuple<std::string, std::string, std::string>> file_triples;
  while (std::getline(input_file, line)) {
    std::stringstream ss(line);
    std::string p_path, q_path, hinge_path;
    getline(ss, p_path, ',');
    getline(ss, q_path, ',');
    getline(ss, hinge_path, ',');
    file_triples.push_back(std::make_tuple(p_path, q_path, hinge_path));
  }
#pragma omp parallel for
  for (const auto &triple : file_triples) {
    std::string p_pdb_id =
        std::get<0>(triple).substr(std::get<0>(triple).find_last_of("/") + 1);
    Eigen::MatrixXd p = openMatrixData(std::get<0>(triple));
    Eigen::MatrixXd q = openMatrixData(std::get<1>(triple));
    std::string hinge_file = std::get<2>(triple);
    std::string hingeIndices = extractHingeIndices(hinge_file);
    int total_residue_length = p.cols();
    if (p.cols() != q.cols()) {
      std::cout << "p length: " << p.cols() << " q length: " << q.cols()
                << std::endl;
      std::cout << "The residue length is different" << std::endl;
      continue;
    }
    if (p.cols() == 0) {
      std::cout << "No data" << std::endl;
      continue;
    }
    std::cout << total_residue_length << std::endl;
    std::vector<double> default_weights;
    for (int i = 0; i < total_residue_length; i++) {
      default_weights.push_back(1.0);
    }
    ConformationPair PQ_pair = MoveToOrigin(p, q, default_weights);
    double rmsd_result = CalcRMSD(PQ_pair.P, PQ_pair.Q, default_weights);
    std::cout << p_pdb_id << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    if (std::string(argv[2]) == "aic") {
      c = 2 * sigma * sigma * (7);
    } else if (std::string(argv[2]) == "bic") {
      c = sigma * sigma * (3 + 4) * log(total_residue_length);
    } else if (std::string(argv[2]) == "lh_aic") {
      c = 2 * sigma * sigma * (7);
    } else if (std::string(argv[2]) == "lh_bic") {
      c = sigma * sigma * (3 + 4) * log(total_residue_length);
    }
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
}
