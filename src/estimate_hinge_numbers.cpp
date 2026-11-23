//go:build ignore

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("unroll-loops")
#endif

#include "rmsdh_new.h"
#include "structure_io.h"

int main(int argc, char *argv[]) try {
  double sigma = 1.55 / sqrt(3);
  double c = 0;
  if (argc != 7) {
    std::cout << "Usage: " << argv[0]
              << " <pdb filename 1> <pdb filename 2> <chain_id 1> <chain id 2> "
                 "<information criterion> <exact or lh>\n";
    return 2;
  }

  std::string pdb_filename1 = argv[1];
  std::string pdb_filename2 = argv[2];
  std::string chain_id1 = argv[3];
  std::string chain_id2 = argv[4];
  const auto valid_chain = [](const std::string &chain) {
    return !chain.empty() && chain.size() <= 8 &&
           std::all_of(chain.begin(), chain.end(), [](unsigned char ch) {
             return (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') ||
                    (ch >= '0' && ch <= '9');
           });
  };
  if (!valid_chain(chain_id1) || !valid_chain(chain_id2) ||
      (std::string(argv[5]) != "aic" && std::string(argv[5]) != "bic") ||
      (std::string(argv[6]) != "exact" && std::string(argv[6]) != "lh")) {
    throw std::invalid_argument("Invalid estimation options");
  }
  const bool is_heuristic = std::string(argv[6]) == "lh";
  const bool pdbNumbering = structure_io::extension(pdb_filename1) != ".cif" ||
                            structure_io::extension(pdb_filename2) != ".cif";
  auto first = structure_io::read(pdb_filename1, chain_id1, pdbNumbering);
  auto second = structure_io::read(pdb_filename2, chain_id2, pdbNumbering);
  std::vector<structure_io::Residue> residues;
  Eigen::Matrix3Xd p(3, std::min(first.size(), second.size())), q(3, p.cols());
  auto a = first.begin(), b = second.begin();
  while (a != first.end() && b != second.end()) {
    if (a->first < b->first) {
      ++a;
    } else if (b->first < a->first) {
      ++b;
    } else {
      p.col(residues.size()) = a->second;
      q.col(residues.size()) = b->second;
      residues.push_back(a->first);
      ++a;
      ++b;
    }
  }
  int total_residue_length = residues.size();
  if (residues.empty())
    throw std::invalid_argument("Structures have no matching CA coordinates");
  p.conservativeResize(3, total_residue_length);
  q.conservativeResize(3, total_residue_length);
  std::cout << "Protein length: " << total_residue_length << std::endl;
  std::vector<double> default_weights(total_residue_length, 1.0);
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
  std::cout << "The number of hinges is: " << hinge_cnt << std::endl;
  std::cout << "Hinge index is: ";
  const auto &indices = rmsdh_hinge_cnt_result.hinge_index_vec;
  for (auto it = indices.rbegin(); it != indices.rend(); ++it) {
    const auto &id = residues.at(*it - 1);
    if (it != indices.rbegin())
      std::cout << " : ";
    std::cout << id.first;
    if (id.second != ' ')
      std::cout << id.second;
  }
  std::cout << std::endl;
  return 0;
} catch (const std::exception &) {
  std::cerr << "Invalid structure or estimation options" << std::endl;
  return 2;
}
