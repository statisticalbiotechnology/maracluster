/******************************************************************************  
  Copyright 2015 Matthew The <matthew.the@scilifelab.se>
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
 ******************************************************************************/
 
#include "MatrixLoader.h"

namespace maracluster {

/* for non-binary input
bool MatrixLoader::initStream(const std::string& matrixFN) {
  matrixStream_.open(matrixFN.c_str());
  if (!matrixStream_.is_open()) {
    std::cerr << "Could not open matrix file " << matrixFN << std::endl;
    return false;
  } else {
    edgesAvailable_ = true;
    return true;
  }
}

// reads in a sparse matrix from a tab delimited file, organised as: fileIdx1 <tab> scannr1 <tab> fileIdx2 <tab> scannr2 <tab> pvalue
bool MatrixLoader::nextEdge(unsigned int& row, unsigned int& col, double& value) {  
  std::string line;
  if (getline(matrixStream_, line)) {
    if (false) {
      std::string filePath1, filePath2;
      std::istringstream iss(line);
      iss >> filePath1 >> row >> filePath2 >> col >> value;
      row = fileList_.getScanId(filePath1, row);
      col = fileList_.getScanId(filePath2, col);
    } else {
      std::istringstream iss(line);
      iss >> row >> col >> value;
    }
    return true;
  } else {
    edgesAvailable_ = false;
    return false;
  }
}
*/

bool MatrixLoader::initStream(const std::string& matrixFN) {
  matrixStream_.open(matrixFN.c_str(), std::ios::in | std::ios::binary);
  if (!matrixStream_.is_open()) {
    std::cerr << "Could not open matrix file " << matrixFN << std::endl;
    return false;
  } else {
    numPvals_ = estimateNumPvals(matrixFN);
    edgesAvailable_ = true;
    return true;
  }
}

// reads in a sparse matrix from a tab delimited file, organised as: fileIdx1 <tab> scannr1 <tab> fileIdx2 <tab> scannr2 <tab> pvalue
bool MatrixLoader::nextEdge(ScanId& row, ScanId& col, double& value) {  
  PvalueTriplet tmp;
  char buffer[sizeof(PvalueTriplet)];
  if (matrixStream_.read(buffer, sizeof(PvalueTriplet))) {
    memcpy(&tmp, buffer, sizeof(tmp));
    row = tmp.scannr1;
    col = tmp.scannr2;
    value = tmp.pval;
    return true;
  } else {
    edgesAvailable_ = false;
    return false;
  }
}

bool MatrixLoader::nextNEdges(unsigned int n, std::vector<PvalueTriplet>& pvec) {  
  PvalueTriplet tmp;
  char buffer[sizeof(PvalueTriplet)];
  unsigned int i = 0;
  while (++i <= n && matrixStream_.read(buffer, sizeof(PvalueTriplet))) {
    memcpy(&tmp, buffer, sizeof(tmp));
    pvec.push_back(tmp);
  }
  
  if (i > n) {
    return true;
  } else {
    edgesAvailable_ = false;
    return false;
  }
}

/*
bool MatrixLoader::initStream(const std::string& matrixFN) {
  numPvals_ = estimateNumPvals(matrixFN);
  mmap_ = boost::iostreams::mapped_file(matrixFN, 
            boost::iostreams::mapped_file::readonly);
  f_ = mmap_.const_data();
  l_ = f_ + mmap_.size();
  
  if (!f_) {
    std::cerr << "Could not open matrix file " << matrixFN << std::endl;
    return false;
  } else {
    edgesAvailable_ = true;
    return true;
  }
}

bool MatrixLoader::nextEdge(unsigned int& row, unsigned int& col, double& value) {  
  PvalueTriplet tmp;
  if (errno == 0 && f_ && f_<(l_-sizeof(tmp)) ) {
    memcpy(&tmp, f_, sizeof(tmp));
    f_ += sizeof(tmp);
    row = tmp.scannr1;
    col = tmp.scannr2;
    value = tmp.pval;
    //std::cerr << row << " " << col << " " << value << std::endl;
    return true;
  } else {
    edgesAvailable_ = false;
    return false;
  }
}
*/

long long MatrixLoader::estimateNumPvals(const std::string& pvalFN) {
  long long fileSize = getFileSize(pvalFN);
  return static_cast<long long>(fileSize / sizeof(PvalueTriplet));
}

long long MatrixLoader::getFileSize(const std::string& pvalFN) {
  std::ifstream in(pvalFN.c_str(), std::ios::ate | std::ios::binary);
  return static_cast<long long>(in.tellg());
}

} /* namespace maracluster */
