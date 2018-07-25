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

#ifndef MARACLUSTER_SPARSEMATRIX_H_
#define MARACLUSTER_SPARSEMATRIX_H_

#include <queue>
#include <vector>
#include <utility>

#include <boost/unordered_map.hpp>

#include "ScanId.h"
#include "SparseEdge.h"

namespace maracluster {

//typedef boost::unordered_map<ScanId, std::map<ScanId, double> > SparseMatrix
class SparseMatrix {
 public:
  SparseMatrix() {
    idxToScanId_.push_back(ScanId(0,0));
    sparseMatrix_.push_back(SparseRow());
  }
  
  inline bool isAlive(const ScanId& si) {
    return isAlive(scanIdToIdx_[si]);
  }
  
  void reserve(const size_t numScans) {
    idxToScanId_.reserve(numScans + 1);
    sparseMatrix_.reserve(numScans + 1);
  }
  
  void insert(const ScanId& s1, const ScanId& s2, const double value) {
    //std::cerr << "Inserting " << s1 << " " << s2 << " " << value << std::endl;
    size_t i1 = getIdx(s1);
    size_t i2 = getIdx(s2);
    insert(i1, i2, value);
  }
  
  void sortRows() {
    for (size_t i1 = 0; i1 < sparseMatrix_.size(); ++i1) {
      std::sort(sparseMatrix_[i1].begin(), sparseMatrix_[i1].end());
    }
  }
  
  // for each element in the SparseRow of s1, check if s2 also has a link to 
  // this element and merge these links if this is the case
  void merge(const ScanId& s1, const ScanId& s2, const ScanId& m, 
             std::priority_queue<SparseEdge>& edgeList) {
    size_t i1 = getIdx(s1);
    size_t i2 = getIdx(s2);
    size_t row = getIdx(m);
    for (size_t idx1 = 0u; idx1 < sparseMatrix_[i1].size(); ++idx1) {
      size_t col = sparseMatrix_[i1][idx1].first;
      
      if (!isAlive(col)) continue;
      
      for (size_t idx2 = 0u; idx2 < sparseMatrix_[i2].size(); ++idx2) {
        if (sparseMatrix_[i2][idx2].first == col) {
#ifdef SINGLE_LINKAGE // TODO: test if this indeed performs single linkage
          double val = std::min(sparseMatrix_[i1][idx1].second, 
                                sparseMatrix_[i2][idx2].second);
#else
          //double value = (sparseMatrix_[i1][idx1].second * n + sparseMatrix_[i2][idx2].second * m)/(n+m); // UPGMA
          double val = std::max(sparseMatrix_[i1][idx1].second, 
                                sparseMatrix_[i2][idx2].second);
#endif
          insert(row, col, val);
          edgeList.push(SparseEdge(val, m, idxToScanId_[col]));
        }
        if (sparseMatrix_[i2][idx2].first >= col) break;
      }
    }
    
    SparseRow empty1, empty2;
    sparseMatrix_[i1].swap(empty1);
    sparseMatrix_[i2].swap(empty2);
  }
  
 protected:
  typedef std::pair<size_t, double> SparseEntry;
  typedef std::vector<SparseEntry> SparseRow;
  std::vector<SparseRow> sparseMatrix_;
  boost::unordered_map<ScanId, size_t> scanIdToIdx_;
  std::vector<ScanId> idxToScanId_;
  
  void insert(size_t i1, size_t i2, const double value) {
    sparseMatrix_[i1].push_back(std::make_pair(i2, value));
    sparseMatrix_[i2].push_back(std::make_pair(i1, value));
  }
  
  inline bool isAlive(const size_t idx) {
    return sparseMatrix_[idx].size() > 0;
  }
  
  size_t getIdx(const ScanId& s1) {
    //std::cerr << "Getting idx for " << s1 << std::endl;
    if (scanIdToIdx_[s1] == 0) {
      scanIdToIdx_[s1] = idxToScanId_.size();
      idxToScanId_.push_back(s1);
      sparseMatrix_.push_back(SparseRow());
    }
    return scanIdToIdx_[s1];
  }
};

} /* namespace maracluster */

#endif /* MARACLUSTER_SPARSEMATRIX_H_ */

