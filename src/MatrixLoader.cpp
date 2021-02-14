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

void MatrixLoader::initVector(const std::vector<PvalueTriplet>& pvec) {
  pvals_ = pvec;
  numPvals_ = pvals_.size();
  edgesAvailable_ = true;
  useFileInput_ = false;
}

bool MatrixLoader::initStream(const std::string& matrixFN) {
  matrixStream_.open(matrixFN.c_str(), std::ios::in | std::ios::binary);
  if (!matrixStream_.is_open()) {
    std::cerr << "Could not open matrix file " << matrixFN << std::endl;
    return false;
  } else {
    numPvals_ = estimateNumPvals(matrixFN);
    edgesAvailable_ = true;
    useFileInput_ = true;
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

void MatrixLoader::nextNEdges(unsigned int n, std::vector<PvalueTriplet>& pvec) {
  if (!useFileInput_) {
    if (offset_ < numPvals_) {
      unsigned int n_min = (std::min)((unsigned int)(numPvals_ - offset_), n);
      pvec.insert(pvec.begin(), pvals_.begin() + offset_, pvals_.begin() + offset_ + n_min);
      offset_ += n_min;
      if (offset_ >= numPvals_) {
        pvals_.clear();
        edgesAvailable_ = false;
      }
    }
  } else {
    PvalueTriplet tmp;
    char buffer[sizeof(PvalueTriplet)];
    unsigned int i = 0;
    while (++i <= n && matrixStream_.read(buffer, sizeof(PvalueTriplet))) {
      memcpy(&tmp, buffer, sizeof(tmp));
      pvec.push_back(tmp);
    }
    
    if (i < n) {
      edgesAvailable_ = false;
    }
  }
}

long long MatrixLoader::estimateNumPvals(const std::string& pvalFN) {
  long long fileSize = getFileSize(pvalFN);
  return static_cast<long long>(fileSize / sizeof(PvalueTriplet));
}

long long MatrixLoader::getFileSize(const std::string& pvalFN) {
  std::ifstream in(pvalFN.c_str(), std::ios::ate | std::ios::binary);
  return static_cast<long long>(in.tellg());
}

} /* namespace maracluster */
