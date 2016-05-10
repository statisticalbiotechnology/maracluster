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
 
#ifndef SPARSE_CLUSTERING_H
#define SPARSE_CLUSTERING_H

#include <map>
#include <utility>
#include <queue>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>

#include <iostream>
#include <fstream>

#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>

#include "MatrixLoader.h"
#include "PvalueTriplet.h"

struct SparseEdge {
  SparseEdge(double _value, ScanId _row, ScanId _col) : 
      value(_value), row(_row), col(_col) {}
  SparseEdge() : value(0.0), row(), col() {}
  
  double value;
  ScanId row, col;
  
  inline bool operator<(const SparseEdge& r2) const {
    return value > r2.value
      || (value == r2.value && row < r2.row)
      || (value == r2.value && row == r2.row && col < r2.col);
  }
  
  bool sameEdge(const SparseEdge& r2) {
    return (row == r2.row && col == r2.col);
  }
};

struct SparseMissingEdge : public SparseEdge {
  SparseMissingEdge(double _value, ScanId _row, ScanId _col, 
             unsigned int _numEdges) : SparseEdge(_value, _row, _col), 
                                       numEdges(_numEdges) {}
  SparseMissingEdge() : SparseEdge(), numEdges(1) {}
  
  unsigned int numEdges;
};

class SparseClustering {
 public:
  SparseClustering() : numTotalEdges_(0), clusterPairFN_(""), 
    edgeLoadingBatchSize_(20000000) /* 20M */, 
    mergeOffset_(3000000000) /* 3G */,
    writeMissingEdges_(false) { }
  
  inline void setClusterPairFN(const std::string& clusterPairFN) { 
    clusterPairFN_ = clusterPairFN;
  }
  
  const std::map<ScanId, std::vector<ScanId> >& getClusters() const { 
    return clusters_; 
  }
  void initMatrix(const std::string& matrixFN);
  
  void doClustering(double cutoff);
  
  void printClusters(std::string& resultFN);
  
  void setMergeOffset(unsigned int mergeOffset) { mergeOffset_ = mergeOffset; }
  
  static bool clusteringUnitTest();
 protected:
  long long numTotalEdges_;
  std::string clusterPairFN_;
  unsigned int edgeLoadingBatchSize_;
  unsigned int mergeOffset_;
  bool writeMissingEdges_;
  
  MatrixLoader matrixLoader_;
  
  std::priority_queue<SparseEdge> edgeList_;
  boost::unordered_map<ScanId, std::map<ScanId, double> > matrix_;
  std::map<ScanId, std::vector<ScanId> > clusters_;
  
  boost::unordered_map<ScanId, ScanId> mergeRoots_;
  boost::unordered_map<ScanId, bool> isAlive_;
  std::vector<SparseMissingEdge> missingEdges_;
  
  void loadNextEdges();
  
  void getClusterMemberships(
    boost::unordered_map<ScanId, ScanId>& clusterMemberships);
  void updateMissingEdges(
    boost::unordered_map<ScanId, ScanId>& clusterMemberships);
  void addNewEdges(
    boost::unordered_map<ScanId, ScanId>& clusterMemberships);
  void pruneEdges();
  
  void addNewEdge(const SparseMissingEdge& edge, size_t& idx);
  void insertEdge(const ScanId row, const ScanId col, const double value);
  void popEdge();
  
  void clusterInit(const ScanId row);
  ScanId getRoot(const ScanId& si);
  
  void joinClusters(const ScanId& minRow, const ScanId& minCol,
    const ScanId& mergeScanId);
  void updateMatrix(const ScanId& minRow, const ScanId& minCol,
    const ScanId& mergeScanId);
  
  void writeMissingEdges(double cutoff);
  
  inline static bool lowerEdge(const SparseMissingEdge& a, 
                               const SparseMissingEdge& b) { 
    return (a.row < b.row) || 
           (a.row == b.row && a.col < b.col) || 
           (a.row == b.row && a.col == b.col && a.value > b.value); 
  }
  
  inline static bool lowerPval(const SparseMissingEdge& a, 
                               const SparseMissingEdge& b) { 
    return (a.value < b.value); 
  }
};

#endif // SPARSE_CLUSTERING_H

