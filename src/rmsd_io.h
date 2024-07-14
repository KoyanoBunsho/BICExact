#ifndef rmsd_io
#define rmsd_io
#include "rmsd_struct.h"

void read_csv(std::vector<std::vector<std::string>> &data,
              std::string csv_path) {
  std::ifstream ifs(csv_path);
  if (!ifs) {
    std::cerr << "Failed to open file." << std::endl;
    return;
  }
  // Skip the column name line
  std::string skipline;
  std::getline(ifs, skipline);

  std::string line;
  while (std::getline(ifs, line)) {
    std::vector<std::string> row;
    size_t pos = 0;
    std::string delimiter = ",";
    while ((pos = line.find(delimiter)) != std::string::npos) {
      std::string token = line.substr(0, pos);
      row.push_back(token);
      line.erase(0, pos + delimiter.length());
    }
    row.push_back(line);
    data.push_back(row);
  }
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

Eigen::MatrixXd openMatrixData(std::string fileToOpen) {
  std::vector<double> matrixEntries;
  std::ifstream matrixDataFile(fileToOpen);
  if (!matrixDataFile) {
    return Eigen::MatrixXd();
  }
  std::string matrixRowString;
  std::string matrixEntry;
  int matrixRowNumber = 0;
  while (getline(matrixDataFile, matrixRowString)) {
    std::stringstream matrixRowStringStream(matrixRowString);
    while (getline(matrixRowStringStream, matrixEntry, ',')) {
      matrixEntries.push_back(stod(matrixEntry));
    }
    matrixRowNumber++;
  }
  if (matrixRowNumber > 0) {
    return Eigen::Map<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        matrixEntries.data(), matrixRowNumber,
        matrixEntries.size() / matrixRowNumber);
  } else {
    return Eigen::MatrixXd();
  }
}

Eigen::MatrixXd openMatrixData2(std::string fileToOpen) {
  std::vector<double> matrixEntries;
  std::ifstream matrixDataFile(fileToOpen);
  if (!matrixDataFile) {
    return Eigen::MatrixXd();
  }
  std::string matrixRowString;
  std::string matrixEntry;
  int matrixRowNumber = 0;
  while (getline(matrixDataFile, matrixRowString)) {
    std::stringstream matrixRowStringStream(matrixRowString);
    while (getline(matrixRowStringStream, matrixEntry, ',')) {
      matrixEntries.push_back(stod(matrixEntry));
    }
    matrixRowNumber++;
  }
  int matrixColumnNumber = 3;

  if (matrixRowNumber > 0) {
    return Eigen::Map<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        matrixEntries.data(), matrixColumnNumber, matrixRowNumber);
  } else {
    return Eigen::MatrixXd();
  }
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
