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

#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <queue>
#include <vector>
#include <utility>

#include <boost/unordered_map.hpp>

#include "ScanId.h"
#include "SparseEdge.h"

//typedef boost::unordered_map<ScanId, std::map<ScanId, double> > SparseMatrix
class SparseMatrix {
 public:
  SparseMatrix() {
    idxToScanId_.push_back(ScanId(0,0));
    sparseMatrix_.push_back(SparseRow());
  }
  
  inline bool isAlive(const ScanId& si) {
    return sparseMatrix_[scanIdToIdx_[si]].size() > 0;
  }
  
  void reserve(const size_t numScans) {
    idxToScanId_.reserve(numScans + 1);
    sparseMatrix_.reserve(numScans + 1);
  }
  
  void insert(const ScanId& s1, const ScanId& s2, const double value) {
    //std::cerr << "Inserting " << s1 << " " << s2 << " " << value << std::endl;
    size_t i1 = getScanIdx(s1);
    size_t i2 = getScanIdx(s2);
    insert(i1, i2, value);
  }
  
  void sortRows() {
    for (size_t i1 = 0; i1 < sparseMatrix_.size(); ++i1) {
      std::sort(sparseMatrix_[i1].begin(), sparseMatrix_[i1].end());
    }
  }
  
  void merge(const ScanId& s1, const ScanId& s2, const ScanId& m, std::priority_queue<SparseEdge>& edgeList) {
    size_t i1 = getScanIdx(s1);
    size_t i2 = getScanIdx(s2);
    size_t row = getScanIdx(m);
    size_t idx2 = 0u;
    for (size_t idx1 = 0u; idx1 < sparseMatrix_[i1].size(); ++idx1) {
      size_t col = sparseMatrix_[i1][idx1].first;
      while (idx2 < sparseMatrix_[i2].size() && sparseMatrix_[i2][idx2].first <= col) {
        if (sparseMatrix_[i2][idx2].first == col && sparseMatrix_[col].size() > 0) {
          double val = std::max(sparseMatrix_[i1][idx1].second, 
                                sparseMatrix_[i2][idx2].second);
          insert(row, col, val);
          edgeList.push(SparseEdge(val, m, idxToScanId_[col]));
        }
        ++idx2;
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
  
  size_t getScanIdx(const ScanId& s1) {
    //std::cerr << "Getting idx for " << s1 << std::endl;
    if (scanIdToIdx_[s1] == 0) {
      scanIdToIdx_[s1] = idxToScanId_.size();
      idxToScanId_.push_back(s1);
      sparseMatrix_.push_back(SparseRow());
    }
    return scanIdToIdx_[s1];
  }
};

#endif

