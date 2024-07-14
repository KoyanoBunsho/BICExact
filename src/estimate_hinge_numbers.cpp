#pragma GCC target("avx2")
#pragma GCC optimize("unroll-loops")

#include "rmsdh_new.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>

class PDBReader {
private:
  std::string pdb_filename;

public:
  PDBReader(const std::string &filename) : pdb_filename(filename) {}

  std::set<int> get_residue_numbers(const std::string &chain_id) const {
    std::ifstream file(pdb_filename);
    std::string line;
    std::set<int> residue_numbers;

    if (file.is_open()) {
      while (getline(file, line)) {
        std::stringstream ss(line);
        std::string field;
        std::getline(ss, field, ' ');
        if (field == "ATOM") {
          std::string atom_name;
          std::string chain;
          int res_num;
          ss >> field;
          ss >> atom_name;
          ss >> field;
          ss >> chain;
          ss >> res_num;
          if (chain == chain_id) {
            residue_numbers.insert(res_num);
          }
        }
      }
      file.close();
    } else {
      std::cout << "Unable to open file: " << pdb_filename << std::endl;
    }

    return residue_numbers;
  }

  std::vector<std::tuple<double, double, double>>
  get_CA_coordinates(const std::set<int> &residue_numbers,
                     const std::string &chain_id) const {
    std::ifstream file(pdb_filename);
    std::string line;
    std::vector<std::tuple<double, double, double>> coordinates;

    if (file.is_open()) {
      while (getline(file, line)) {
        std::stringstream ss(line);
        std::string field;
        std::getline(ss, field, ' ');
        if (field == "ATOM") {
          std::string atom_name;
          std::string chain;
          int res_num;
          double x, y, z;
          ss >> field;
          ss >> atom_name;
          ss >> field;
          ss >> chain;
          ss >> res_num;
          ss >> x >> y >> z;
          if (atom_name == "CA" && chain == chain_id &&
              residue_numbers.find(res_num) != residue_numbers.end()) {
            coordinates.push_back(std::make_tuple(x, y, z));
          }
        }
      }
      file.close();
    } else {
      std::cout << "Unable to open file: " << pdb_filename << std::endl;
    }

    return coordinates;
  }
};

std::set<int> intersect_residue_numbers(const std::set<int> &set1,
                                        const std::set<int> &set2) {
  std::set<int> intersection;
  std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
                        std::inserter(intersection, intersection.begin()));
  return intersection;
}

Eigen::MatrixXd
convert_to_matrix(const std::vector<std::tuple<double, double, double>> &vec) {
  Eigen::MatrixXd matrix(vec.size(), 3);
  for (int i = 0; i < (int)vec.size(); ++i) {
    matrix(i, 0) = std::get<0>(vec[i]);
    matrix(i, 1) = std::get<1>(vec[i]);
    matrix(i, 2) = std::get<2>(vec[i]);
  }
  return matrix;
}

int main(int argc, char *argv[]) {
  double sigma = 1.55 / sqrt(3);
  double c = 0;
  bool is_heuristic = false;
  if (argc != 7) {
    std::cout << "Usage: " << argv[0]
              << " <pdb filename 1> <pdb filename 2> <chain_id 1> <chain id 2> "
                 "<information criterion> <exact or lh>\n";
    return 1;
  }

  std::string pdb_filename1 = argv[1];
  std::string pdb_filename2 = argv[2];
  std::string chain_id1 = argv[3];
  std::string chain_id2 = argv[4];
  if (std::string(argv[6]) == "exact") {
    is_heuristic = false;
  } else if (std::string(argv[6]) == "lh") {
    is_heuristic = true;
  }
  PDBReader reader1(pdb_filename1);
  PDBReader reader2(pdb_filename2);
  std::set<int> res_numbers1 = reader1.get_residue_numbers(chain_id1);
  std::set<int> res_numbers2 = reader2.get_residue_numbers(chain_id2);
  std::set<int> intersected_res_numbers =
      intersect_residue_numbers(res_numbers1, res_numbers2);
  std::vector<std::tuple<double, double, double>> coordinatesP =
      reader1.get_CA_coordinates(intersected_res_numbers, chain_id1);
  std::vector<std::tuple<double, double, double>> coordinatesQ =
      reader2.get_CA_coordinates(intersected_res_numbers, chain_id2);
  Eigen::MatrixXd p, q;
  p = convert_to_matrix(coordinatesP).transpose();
  q = convert_to_matrix(coordinatesQ).transpose();
  int total_residue_length = p.cols();
  std::cout << "Protein length: " << total_residue_length << std::endl;
  std::vector<double> default_weights;
  for (int i = 0; i < total_residue_length; i++) {
    default_weights.push_back(1.0);
  }
  ConformationPair PQ_pair = MoveToOrigin(p, q, default_weights);
  double rmsd_result = CalcRMSD(PQ_pair.P, PQ_pair.Q, default_weights);
  std::cout << "RMSD: " << rmsd_result << std::endl;
  if (std::string(argv[5]) == "aic") {
    c = 2 * sigma * sigma * (7);
  } else if (std::string(argv[5]) == "bic") {
    c = sigma * sigma * (3 + 4) * log(total_residue_length);
  }
  ProteinRMSDhinge rmsdh_calculator(PQ_pair.P, PQ_pair.Q, c);
  RMSDhHingeCnt rmsdh_hinge_cnt_result;
  if (is_heuristic) {
    rmsdh_hinge_cnt_result = rmsdh_calculator.CalcLinearRMSDh();
  } else {
    rmsdh_hinge_cnt_result = rmsdh_calculator.CalcRMSDh();
  }
  double rmsdh_result = rmsdh_hinge_cnt_result.rmsdh_result;
  std::cout << "RMSDhk is: " << rmsdh_result << std::endl;
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
  std::cout << "The number of hinges is: " << hinge_cnt << std::endl;
  std::cout << "Hinge index is: " << hinge_index << std::endl;
  return 0;
}
