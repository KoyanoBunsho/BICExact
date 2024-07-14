#pragma GCC target("avx2")
#pragma GCC optimize("unroll-loops")

#include "rmsd_io.h"
#include "rmsd_struct.h"
#include "rmsdh_new.h"

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << "Input parameter" << std::endl;
    return 0;
  }
  double sigma = std::stod(argv[1]) / sqrt(3);
  bool is_heuristic = false;
  std::string coord_path;
  std::string database_path;
  std::string save_database_name;
  std::string save_method_name;
  double c = 0;
  if (std::string(argv[2]) == "dyndom") {
    coord_path = "coord_csv_dyndom/";
    database_path = "DynDon_database.csv";
    save_database_name = "dyndom";
  } else if (std::string(argv[2]) == "shibuya") {
    coord_path = "coord_csv/";
    database_path = "pdb_12_with_hinges.csv";
    save_database_name = "shibuya";
  } else {
    std::cerr << "Choose dyndom or shibuya" << std::endl;
  }
  if (std::string(argv[3]) == "aic") {
    save_method_name = "rmsdh_aic";
  } else if (std::string(argv[3]) == "bic") {
    save_method_name = "rmsdh_bic";
  } else if (std::string(argv[3]) == "lh_aic") {
    save_method_name = "fast_rmsdh_aic";
    is_heuristic = true;
  } else if (std::string(argv[3]) == "lh_bic") {
    save_method_name = "fast_rmsdh_bic";
    is_heuristic = true;
  }
  std::ofstream myfile;
  std::vector<std::vector<std::string>> pdb_chain_data;
  read_csv(pdb_chain_data, database_path);
  std::vector<std::pair<std::string, std::string>> pdb_pair_vec;
  for (const auto &pdb_chain : pdb_chain_data) {
    pdb_pair_vec.push_back(std::make_pair(pdb_chain[0], pdb_chain[1]));
  }
  std::string save_name = "rmsdh_result/" + save_database_name + "_" +
                          save_method_name + "_" + std::to_string(sigma) +
                          ".csv";
  myfile.open(save_name);
  myfile << "p_pdb_id,q_pdb_id,Residue "
            "length,RMSD,RMSDh,c,hinge_cnt,hinge_index,exec_time (s)"
         << std::endl;
  for (int i = 0; i < (int)pdb_pair_vec.size(); i++) {
    std::string p_pdb_id, q_pdb_id;
    std::string p_pdb_chain_id = pdb_pair_vec[i].first;
    std::string q_pdb_chain_id = pdb_pair_vec[i].second;
    std::transform(p_pdb_chain_id.begin(), p_pdb_chain_id.begin() + 4,
                   std::back_inserter(p_pdb_id), ::tolower);
    std::string p_chain_id = p_pdb_chain_id.substr(5, 1);
    std::transform(q_pdb_chain_id.begin(), q_pdb_chain_id.begin() + 4,
                   std::back_inserter(q_pdb_id), ::tolower);
    std::string q_chain_id = q_pdb_chain_id.substr(5, 1);
    Eigen::MatrixXd p, q;
    p = openMatrixData(coord_path + "coord_" + p_pdb_id + "_" + p_chain_id +
                       "_" + q_pdb_id + "_" + q_chain_id + ".csv");
    std::cout << coord_path + "coord_" + p_pdb_id + "_" + p_chain_id + "_" +
                     q_pdb_id + "_" + q_chain_id + ".csv"
              << std::endl;
    q = openMatrixData(coord_path + "coord_" + q_pdb_id + "_" + q_chain_id +
                       "_" + p_pdb_id + "_" + p_chain_id + ".csv");
    int total_residue_length = p.cols();
    if (p.cols() != q.cols()) {
      std::cout << "p length: " << p.cols() << " q length: " << q.cols()
                << std::endl;
      std::cout << "The residue length is different" << std::endl;
      continue;
    }
    std::cout << total_residue_length << std::endl;
    if (total_residue_length == 0) {
      std::cout << "No data" << std::endl;
      continue;
    }
    std::vector<double> default_weights;
    for (int i = 0; i < total_residue_length; i++) {
      default_weights.push_back(1.0);
    }
    ConformationPair PQ_pair = MoveToOrigin(p, q, default_weights);
    double rmsd_result = CalcRMSD(PQ_pair.P, PQ_pair.Q, default_weights);
    myfile << p_pdb_chain_id << "," << q_pdb_chain_id << ",";
    myfile << total_residue_length << ",";
    myfile << rmsd_result << ",";
    std::cout << p_pdb_id << ", " << q_pdb_id << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    if (std::string(argv[3]) == "aic") {
      c = 2 * sigma * sigma * (7);
    } else if (std::string(argv[3]) == "bic") {
      c = sigma * sigma * (3 + 4) * log(total_residue_length);
    } else if (std::string(argv[3]) == "lh_aic") {
      c = 2 * sigma * sigma * (7);
    } else if (std::string(argv[3]) == "lh_bic") {
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
    myfile << rmsdh_result << "," << c << "," << hinge_cnt << "," << hinge_index
           << "," << exec_time_s << std::endl;
  }
  myfile.close();
}
