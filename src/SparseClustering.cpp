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
 
#include "SparseClustering.h"

void SparseClustering::initMatrix(const std::string& matrixFN) {
  matrixLoader_.initStream(matrixFN);
}

void SparseClustering::clusterInit(const ScanId& row) {
#pragma omp critical (new_cluster)
  {
    isAlive_[row] = true;
    if (clusters_[row].size() == 0) {
      clusters_[row].push_back(row);
    }
  }
}

void SparseClustering::loadNextEdges() {
  boost::unordered_map<ScanId, ScanId> clusterMemberships;
  getClusterMemberships(clusterMemberships);
#ifndef SINGLE_LINKAGE
  updateMissingEdges(clusterMemberships);
#endif
  addNewEdges(clusterMemberships);
  
  size_t beforeSize = edgeList_.size();

#ifndef SINGLE_LINKAGE
  pruneEdges();
#endif
  
  matrix_.sortRows();
  
  std::cerr << "  Loaded new edges: new: " << edgeList_.size() - beforeSize << 
      ", total: " << numTotalEdges_ << "/" << matrixLoader_.numPvals_ << 
      " (" << numTotalEdges_*100/matrixLoader_.numPvals_ << "%)." << std::endl;
}

void SparseClustering::getClusterMemberships(
    boost::unordered_map<ScanId, ScanId>& clusterMemberships) {
  std::cerr << "  Creating cluster membership map." << std::endl;
  std::map<ScanId, std::vector<ScanId> >::const_iterator it;
  for (it = clusters_.begin(); it != clusters_.end(); ++it) {
    if (it->second.size() > 0) {
      for (std::vector<ScanId>::const_iterator it2 = it->second.begin(); 
              it2 != it->second.end(); ++it2) {
        clusterMemberships[*it2] = it->first;
      }
    }
  }
}

void SparseClustering::updateMissingEdges(
    boost::unordered_map<ScanId, ScanId>& clusterMemberships) {
  std::cerr << "  Updating incomplete edges (" << missingEdges_.size() << ")." << std::endl;
#pragma omp parallel for schedule(dynamic, 10000) 
  for (int i = 0; i < missingEdges_.size(); ++i) {
    ScanId row = missingEdges_[i].row;
    ScanId col = missingEdges_[i].col;
    ScanId r = getRoot(row);
    ScanId c = getRoot(col);
    boost::unordered_map<ScanId, ScanId>::iterator it;
    if ((it = clusterMemberships.find(r)) != clusterMemberships.end()) {
      r = it->second;
    }
    if ((it = clusterMemberships.find(c)) != clusterMemberships.end()) {
      c = it->second;
    }
    if (r != row || c != col) {
      ScanId s1 = (std::min)(r,c);
      ScanId s2 = (std::max)(r,c);
      
      missingEdges_[i].row = s1;
      missingEdges_[i].col = s2;
    }
  }
}

void SparseClustering::loadEdges(std::vector<PvalueTriplet>& pvec) {
  pvec.reserve(edgeLoadingBatchSize_);
  matrixLoader_.nextNEdges(edgeLoadingBatchSize_, pvec);
}

void SparseClustering::addNewEdges(
    boost::unordered_map<ScanId, ScanId>& clusterMemberships) {
  std::cerr << "  Loading new edges." << std::endl;
  
  std::vector<PvalueTriplet> pvec;
  loadEdges(pvec);
  
  size_t numNewEdges = pvec.size();
  size_t insertOffset = missingEdges_.size();
  numTotalEdges_ += numNewEdges;
#ifndef SINGLE_LINKAGE
  missingEdges_.resize(missingEdges_.size() + numNewEdges);
#endif
#pragma omp parallel for schedule(dynamic, 10000) 
  for (int i = 0; i < numNewEdges; ++i) {
    ScanId row = pvec[i].scannr1;
    ScanId col = pvec[i].scannr2;
    
    boost::unordered_map<ScanId, ScanId>::iterator it;
    if ((it = clusterMemberships.find(row)) != clusterMemberships.end()) {
      row = it->second;
    } else {
      clusterInit(row);
    }
    
    if ((it = clusterMemberships.find(col)) != clusterMemberships.end()) {
      col = it->second;
    } else {
      clusterInit(col);
    }
    
    ScanId s1 = (std::min)(row, col);
    ScanId s2 = (std::max)(row, col);
#ifdef SINGLE_LINKAGE
    if (row != col) {
      size_t idx = 0;
      addNewEdge(SparseMissingEdge(pvec[i].pval, s1, s2, clusters_[s1].size()*clusters_[s2].size()), idx);
    }
#else
    missingEdges_[insertOffset + i] = SparseMissingEdge(pvec[i].pval, s1, s2, 1u);
#endif
  }
}

// Sorts the missing edges and collapses duplicate edges by summing the number 
// of edges. If the number of edges reaches complete linkage, the edge is added 
// to the priority queue
void SparseClustering::pruneEdges() {
  std::cerr << "  Sorting " << missingEdges_.size() << " edges." << std::endl;
  std::sort(missingEdges_.begin(), missingEdges_.end(), lowerEdge);
  
  std::cerr << "  Adding new edges." << std::endl;
  size_t batchSize = 100000;
  size_t collapsedIdx = 0;
#pragma omp parallel for schedule(dynamic, 1) ordered 
  for (int i = 0; i < missingEdges_.size(); i += batchSize) {    
    SparseMissingEdge lastEdge = missingEdges_[i];
    lastEdge.numEdges = 0;
    size_t insertionIdx = i;
    size_t startIdx = i;
    for (size_t j = i; j < std::min(i + batchSize, missingEdges_.size()); ++j) {
      if (lastEdge.sameEdge(missingEdges_[j])) {
        lastEdge.numEdges += missingEdges_[j].numEdges;
      } else {
        addNewEdge(lastEdge, insertionIdx);
        lastEdge = missingEdges_[j];
      }
    }
    addNewEdge(lastEdge, insertionIdx);
  #pragma omp ordered
    {
      // make sure that identical edges spilling over the batches are joined as well
      if (collapsedIdx > 0 && missingEdges_[startIdx].sameEdge(missingEdges_[collapsedIdx-1])) {
        missingEdges_[--collapsedIdx].numEdges += missingEdges_[startIdx++].numEdges;
        addNewEdge(missingEdges_[collapsedIdx], collapsedIdx);
      }
      for (size_t j = startIdx; j < insertionIdx; ++j) {
        missingEdges_[collapsedIdx++] = missingEdges_[j]; // swapping is not faster
      }
    }
  }
  
  missingEdges_.resize(collapsedIdx);
}

void SparseClustering::addNewEdge(const SparseMissingEdge& edge, size_t& idx) {
  if (edge.row != edge.col) {
    if (edge.numEdges == clusters_[edge.row].size()*clusters_[edge.col].size()) {
      insertEdge(edge.row, edge.col, edge.value);
    }
#ifndef SINGLE_LINKAGE
    missingEdges_[idx++] = edge;
#endif
  }
}

void SparseClustering::insertEdge(const ScanId& row, const ScanId& col, 
                                  const double value) {
#pragma omp critical (pq_edges_insert)
  {
    edgeList_.push(SparseEdge(value, row, col));
    matrix_.insert(row, col, value);
  }
}

bool SparseClustering::edgesLeft() {
  return matrixLoader_.hasEdgesAvailable();
}

void SparseClustering::popEdge() {
  edgeList_.pop();
  if (edgeList_.size() == 0 && edgesLeft()) {
    loadNextEdges();
  }
}

void SparseClustering::joinClusters(const ScanId& minRow, const ScanId& minCol,
    const ScanId& mergeScanId) { 
  clusters_[mergeScanId].swap(clusters_[minRow]);
  clusters_[mergeScanId].insert(clusters_[mergeScanId].end(), 
      clusters_[minCol].begin(), clusters_[minCol].end());
  
  clusters_.erase(minRow);
  clusters_.erase(minCol);
}

void SparseClustering::updateMatrix(const ScanId& minRow, const ScanId& minCol,
    const ScanId& mergeScanId) {
  matrix_.merge(minRow, minCol, mergeScanId, edgeList_);
  /*
  typedef std::pair<ScanId, double> ColValPair;
  BOOST_FOREACH (const ColValPair& colValPair, matrix_.getRow(minRow)) {
    ScanId col = colValPair.first;
    if (col == minCol || col == minRow || !isAlive_[col]) continue;
    if (matrix_[minCol].find(col) != matrix_[minCol].end()) {
      //double value = (colValPair.second * n + matrix_[minCol][col] * m)/(n+m); // UPGMA
#ifdef SINGLE_LINKAGE
      double value = (std::min)(colValPair.second, matrix_[minCol][col]); // single linkage
#else
      double value = (std::max)(colValPair.second, matrix_[minCol][col]); // complete linkage
#endif
      insertEdge(mergeScanId, col, value);
    }
  }
  
  matrix_.erase(minRow);
  matrix_.erase(minCol);
  */
}

// Based on http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2718652/
void SparseClustering::doClustering(double cutoff) {  
  std::cerr << "Starting MinHeap clustering" << std::endl;
  
  // decide if we write the results to a file or stdout
  std::ofstream resultFNStream;
  bool writeTree = false;
  if (clusterPairFN_.size() > 0) {
    resultFNStream.open(clusterPairFN_.c_str());
    writeTree = true;
  }
  
  time_t startTime, elapsedTime;
  time(&startTime);
  clock_t startClock = clock(), elapsedClock;
  
  if (edgesLeft()) {
    loadNextEdges();
  } else {
    std::cerr << "Could not read edges from input file." << std::endl;
    return;
  }
  
  unsigned int mergeCnt = 0u;
  while (!edgeList_.empty() && edgeList_.top().value < cutoff) {
    SparseEdge minEdge = edgeList_.top();
    if (isAlive_[minEdge.row] && isAlive_[minEdge.col]) {
      isAlive_[minEdge.row] = false;
      isAlive_[minEdge.col] = false;
      
      if (mergeCnt % 10000 == 0) {
        std::cerr << "It. " << mergeCnt << ": minRow = " << minEdge.row 
                  << ", minCol = " << minEdge.col 
                  << ", minEl = " << minEdge.value << std::endl;
      }
      
      ScanId minRowRoot = getRoot(minEdge.row);
      ScanId mergeScanId(mergeOffset_, mergeCnt++);
      setRoot(mergeScanId, minRowRoot);
      isAlive_[mergeScanId] = true;
      
      if (writeTree) {
        ScanId minColRoot = getRoot(minEdge.col);
        PvalueTriplet tmp(minRowRoot, minColRoot, minEdge.value);
        resultFNStream << tmp << "\n";
      }
      
      joinClusters(minEdge.row, minEdge.col, mergeScanId);
      updateMatrix(minEdge.row, minEdge.col, mergeScanId);
    }
    popEdge();
  }
  
  std::cerr << "Finished MinHeap clustering" << std::endl;
  
  if (writeMissingEdges_) writeMissingEdges(cutoff);
  
  time(&elapsedTime);
  double diff = difftime(elapsedTime, startTime);
  
  unsigned int timeElapsedMin = static_cast<unsigned int>(diff/60);
  unsigned int timeElapsedSecMod = 
      static_cast<unsigned int>(diff - timeElapsedMin * 60);
  
  elapsedClock = clock();
  double elapsedTimeSec = (elapsedClock - startClock) / (double)CLOCKS_PER_SEC;
  std::cerr << "  Elapsed time: " << elapsedTimeSec << " cpu seconds " <<
               "or " << timeElapsedMin << " min " << timeElapsedSecMod << 
               " sec wall time." << std::endl;
}

void SparseClustering::writeMissingEdges(double cutoff) {
  std::cerr << "Writing missing edges." << std::endl;
  loadNextEdges();
  std::sort(missingEdges_.begin(), missingEdges_.end(), lowerPval);
  std::string clusterPairMissingFN = clusterPairFN_ + ".missing_edges.tsv";
  std::ofstream resultMissingFNStream(clusterPairMissingFN.c_str());
  if (resultMissingFNStream.is_open()) {
    BOOST_FOREACH (SparseMissingEdge& sme, missingEdges_) {
      if (sme.value < cutoff - 5.0) {
        sme.row = getRoot(sme.row);
        sme.col = getRoot(sme.col);
        resultMissingFNStream << sme.row << "\t" << sme.col << "\t" << sme.value << "\t" << sme.numEdges << "\n";
      } else {
        break;
      }
    }
  }
}

void SparseClustering::printClusters(std::string& resultFN) {
  // write clustering output to file or stdout
  std::ofstream resultFNStream;
  if (resultFN.size() > 0) {
    resultFNStream.open(resultFN.c_str());
  }
  std::ostream& resultStream = (resultFN.size() > 0) ? resultFNStream : std::cout;
  
  for (std::map<ScanId, std::vector<ScanId> >::const_iterator it = clusters_.begin(); it != clusters_.end(); ++it) {
    if (it->second.size() > 0) {
      for (std::vector<ScanId>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
        resultStream << *it2 << "\t";
      }
      resultStream << std::endl;
    }
  }
}

bool SparseClustering::clusteringUnitTest() {
  
  SparseClustering matrix;
  
  ScanId s1(0,1), s2(0,2), s3(0,3), s4(0,4);
  matrix.clusterInit(s1);
  matrix.clusterInit(s2);
  matrix.clusterInit(s3);
  matrix.clusterInit(s4);
  
  size_t idx = 0u;
  matrix.addNewEdge(SparseMissingEdge(-20.0, s1, s1, 1), idx);
  matrix.addNewEdge(SparseMissingEdge(-2.0, s1, s2, 1), idx);
  matrix.addNewEdge(SparseMissingEdge(-6.0, s1, s3, 1), idx);
  matrix.addNewEdge(SparseMissingEdge(-6.0, s1, s4, 1), idx);
  matrix.addNewEdge(SparseMissingEdge(-1.0, s2, s3, 1), idx);
  matrix.addNewEdge(SparseMissingEdge(-10.0, s2, s4, 1), idx);
  matrix.addNewEdge(SparseMissingEdge(-4.0, s3, s4, 1), idx);
  
  matrix.doClustering(-0.9);
  
  const std::map<ScanId, std::vector<ScanId> >& clusters_ = matrix.getClusters();
  
  /*
  typedef std::pair<ScanId, std::vector<ScanId> > ClusterRow;
  BOOST_FOREACH(ClusterRow row, clusters_) {
    std::cerr << row.first;
    BOOST_FOREACH(ScanId col, row.second) {
      std::cerr << " " << col;
    }
    std::cerr << std::endl;
  }
  */
  
  
  /* UPGMA
  Iteration 0: minRow = 2, minCol = 4, minEl = -10
  Iteration 1: minRow = 1, minCol = 3, minEl = -6
  Iteration 2: minRow = 1, minCol = 2, minEl = -3.25
  */
  
  std::vector<ScanId> clusterRow = clusters_.begin()->second;
  if (clusterRow[0] == s1 && clusterRow[1] == s3 
      && clusterRow[2] == s2 && clusterRow[3] == s4) {
    return true;
  } else {
    std::cerr << "Nr elements = " << clusterRow.size() << ": ";
    BOOST_FOREACH (ScanId scannr, clusterRow) {
      std::cerr << " " << scannr;
    }
    std::cerr << std::endl;
    return false;
  }
}


