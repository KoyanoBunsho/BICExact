#ifndef rmsd_struct
#define rmsd_struct
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Geometry"
#include <cmath>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct RMSDhHingeCnt {
  double rmsdh_result;
  int hinge_cnt;
  std::vector<int> hinge_index_vec;
};
struct AblationResult {
  int iter_num;
  double rmsdh_result;
  int hinge_cnt;
  std::vector<int> hinge_index_vec;
};

struct ConformationPair {
  Eigen::Matrix3Xd P, Q;
};

struct AcNumberPQpdb {
  std::string ac_number;
  std::string pdb_id;
};
#endif
