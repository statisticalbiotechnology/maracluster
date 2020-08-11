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
 
#ifndef MARACLUSTER_SPARSECLUSTERING_H_
#define MARACLUSTER_SPARSECLUSTERING_H_

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
#include <boost/unordered/unordered_map.hpp>

#include "MatrixLoader.h"
#include "PvalueTriplet.h"
#include "SparseEdge.h"
#include "SparseMatrix.h"

namespace maracluster {

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
  
  virtual void doClustering(double cutoff);
  
  void printClusters(std::string& resultFN);
  
  void setMergeOffset(unsigned int mergeOffset) { mergeOffset_ = mergeOffset; }
  
  void reserve(const size_t numScans) {
    matrix_.reserve(numScans + 1);
  }
  
  static bool clusteringUnitTest();
 protected:
  long long numTotalEdges_;
  std::string clusterPairFN_;
  unsigned int edgeLoadingBatchSize_;
  unsigned int mergeOffset_;
  bool writeMissingEdges_;
  
  MatrixLoader matrixLoader_;
  
  std::priority_queue<SparseEdge> edgeList_;
  SparseMatrix matrix_;
  std::map<ScanId, std::vector<ScanId> > clusters_;
  
  boost::unordered_map<ScanId, ScanId> mergeRoots_;
  std::vector<SparseMissingEdge> missingEdges_;
  
  virtual void loadNextEdges();
  void loadEdges(std::vector<PvalueTriplet>& pvec);
  virtual bool edgesLeft();
  
  virtual bool verifyEdgeBoth(SparseEdge& minEdge) { return true; }
  virtual bool verifyEdgeSingle(SparseEdge& minEdge) { return true; }
  
  void getClusterMemberships(
    boost::unordered_map<ScanId, ScanId>& clusterMemberships);
  void updateMissingEdges(
    boost::unordered_map<ScanId, ScanId>& clusterMemberships);
  void addNewEdges(
    boost::unordered_map<ScanId, ScanId>& clusterMemberships);
  void pruneEdges();
  
  void addNewEdge(const SparseMissingEdge& edge, size_t& idx);
  void insertEdge(const ScanId& row, const ScanId& col, const double value);
  void popEdge();
  
  void clusterInit(const ScanId& row);
  
  ScanId getRoot(const ScanId& si) {
    if (si.fileIdx >= mergeOffset_) {
      return mergeRoots_[si];
    } else {
      return si;
    } 
  }
  void setRoot(const ScanId& si, const ScanId& root) { mergeRoots_[si] = root; }
  
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

} /* namespace maracluster */

#endif /* MARACLUSTER_SPARSECLUSTERING_H_ */

