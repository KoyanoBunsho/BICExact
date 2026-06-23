#pragma GCC target("avx2")
#pragma GCC optimize("unroll-loops")

#include "rmsd_io.h"
#include "rmsd_struct.h"
#include "rmsdh_new.h"
#include <experimental/filesystem>
#include <regex>
#include <sstream>

namespace fs = std::experimental::filesystem;

int main(int argc, char **argv) {
  bool is_heuristic = false;
  std::string save_method_name;
  std::string coord_path = "coord_csv_simulation/";
  std::ofstream myfile;
  std::string simulation_data_path = "simulation_data/";
  std::string simulation_data_info_path = "simulation_data_info/";
  std::string input_file_path = std::string(argv[1]);
  std::ifstream input_file(input_file_path);
  std::string save_name = "rmsdh_result/simulation_delta_g_" + input_file_path;
  myfile.open(save_name);
  myfile << "p_pdb_id,Residue length,exec_time (s)" << std::endl;
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
  double c = 100;
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
    ProteinRMSDhinge rmsdh_calculator(PQ_pair.P, PQ_pair.Q, c);
    double delta_g = rmsdh_calculator.CalcDeltaG();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> exec_time_ms = end - start;
    double exec_time_s = exec_time_ms.count() / 1000.0;
    std::cout << exec_time_s << " s" << std::endl;
    myfile << p_pdb_id << ",";
    myfile << total_residue_length << ",";
    myfile << delta_g << "," << exec_time_s << std::endl;
  }
  myfile.close();
}
