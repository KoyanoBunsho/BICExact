#ifndef rmsd_io
#define rmsd_io
#include "rmsd_struct.h"

#include <algorithm>
#include <cctype>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <tuple>

inline std::vector<std::vector<std::string>> read_csv(const std::string &path) {
  std::ifstream input(path);
  std::string line;
  if (!std::getline(input, line))
    throw std::invalid_argument("Cannot read CSV header");
  std::vector<std::vector<std::string>> rows;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    if (line.empty())
      continue;
    std::vector<std::string> row;
    std::size_t start = 0, end;
    do {
      end = line.find(',', start);
      row.push_back(line.substr(start, end - start));
      start = end + 1;
    } while (end != std::string::npos);
    rows.push_back(std::move(row));
  }
  if (input.bad() || rows.empty())
    throw std::invalid_argument("Empty or unreadable CSV");
  return rows;
}

inline std::vector<std::pair<std::string, std::string>>
read_structure_pairs(const std::string &path) {
  std::vector<std::pair<std::string, std::string>> pairs;
  auto valid = [](const std::string &id) {
    return id.size() == 6 && id[4] == '_' &&
           std::all_of(id.begin(), id.begin() + 4,
                       [](unsigned char c) { return std::isalnum(c); }) &&
           std::isalnum(static_cast<unsigned char>(id[5]));
  };
  for (const auto &row : read_csv(path)) {
    if (row.size() < 2 || !valid(row[0]) || !valid(row[1]))
      throw std::invalid_argument("Invalid PDB and chain ID in dataset");
    pairs.emplace_back(row[0], row[1]);
  }
  return pairs;
}

inline std::vector<std::tuple<std::string, std::string, std::string>>
read_simulation_files(const std::string &path) {
  std::vector<std::tuple<std::string, std::string, std::string>> files;
  for (const auto &row : read_csv(path)) {
    if (row.size() < 3 || row[0].empty() || row[1].empty() || row[2].empty())
      throw std::invalid_argument("Expected three paths per simulation record");
    files.emplace_back(row[0], row[1], row[2]);
  }
  return files;
}

inline double read_sigma(const std::string &value) {
  std::size_t end;
  const double sigma = std::stod(value, &end) / std::sqrt(3.0);
  if (end != value.size() || sigma < 0 ||
      !std::isfinite(7 * sigma * sigma *
                     std::log(std::numeric_limits<int>::max())))
    throw std::invalid_argument("Invalid sigma");
  return sigma;
}

template <typename T>
bool getFileContent(std::string fileName, std::vector<T> &Flexibility) {
  std::ifstream in(fileName.c_str());
  if (!in) {
    std::cerr << "Cannot open the File : " << fileName << std::endl;
    return false;
  }
  T val;
  while (in >> val) {
    Flexibility.push_back(val);
  }
  in.close();
  return true;
}

inline Eigen::MatrixXd read_coordinate_matrix(const std::string &path,
                                              bool triples) {
  std::ifstream file(path);
  std::vector<double> values;
  std::string line;
  std::size_t rows = 0, columns = 0;
  while (std::getline(file, line)) {
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    if ((!triples && rows == 3) || line.empty() || line.back() == ',')
      return {};
    std::istringstream row(line);
    std::string field;
    std::size_t count = 0;
    while (std::getline(row, field, ',')) {
      std::size_t end;
      double value;
      try {
        value = std::stod(field, &end);
      } catch (const std::invalid_argument &) {
        return {};
      } catch (const std::out_of_range &) {
        return {};
      }
      if (field.find_first_not_of(" \t", end) != std::string::npos ||
          !std::isfinite(value) || std::abs(value) > 1e6)
        return {};
      values.push_back(value);
      ++count;
    }
    if (rows == 0)
      columns = count;
    if (count != columns || (triples && count != 3))
      return {};
    ++rows;
  }
  if (file.bad() || values.empty() || (!triples && rows != 3))
    return {};
  using RowMatrix =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  return Eigen::Map<const RowMatrix>(values.data(), 3, values.size() / 3);
}

inline Eigen::MatrixXd openMatrixData(const std::string &path) {
  return read_coordinate_matrix(path, false);
}

inline Eigen::MatrixXd openMatrixData2(const std::string &path) {
  return read_coordinate_matrix(path, true);
}

bool getFileStrContent(
    std::string fileName,
    std::vector<std::pair<std::string, std::string>> &pq_pair) {
  std::ifstream in(fileName.c_str());
  if (!in) {
    std::cerr << "Cannot open the File : " << fileName << std::endl;
    return false;
  }
  std::string val_p, val_q;
  while (in >> val_p >> val_q) {
    pq_pair.push_back(make_pair(val_p, val_q));
  }
  in.close();
  return true;
}
#endif
