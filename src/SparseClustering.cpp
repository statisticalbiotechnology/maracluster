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

void SparseClustering::minInsert(const ScanId row, const ScanId col, 
                                 const double value) {
  //edgeList_.push(SparseEdge(value, row, col));
  //inMatrix_[row][col] = (std::min)(value, inMatrix_[row][col]);
  //matrix[col][row] = (std::min)(value, matrix[col][row]);
  matrix_[row][col] = value;
  matrix_[col][row] = value;
  // initialize the cluster scannr list if not initialized previously
  clusterInit(row, col);
}

void SparseClustering::clusterInit(const ScanId row, const ScanId col) {
  clusterInit(row);
  clusterInit(col);
}

void SparseClustering::clusterInit(const ScanId row) {
#pragma omp critical (new_cluster)
  {
    isAlive_[row] = true;
    if (clusters_[row].size() == 0) {
      clusters_[row].push_back(row);
    }
  }
}

void SparseClustering::loadNextEdges(
    boost::unordered_map<ScanId, ScanId>& mergeRoots) {
  boost::unordered_map<ScanId, ScanId> clusterMemberships;
  getClusterMemberships(clusterMemberships);
#ifndef SINGLE_LINKAGE
  updateMissingEdges(clusterMemberships, mergeRoots);
#endif
  addNewEdges(clusterMemberships);
  
  size_t beforeSize = edgeList_.size();

#ifndef SINGLE_LINKAGE
  pruneEdges();
#endif
  
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
    boost::unordered_map<ScanId, ScanId>& clusterMemberships, 
    boost::unordered_map<ScanId, ScanId>& mergeRoots) {
  std::cerr << "  Updating incomplete edges (" << missingEdges_.size() << ")." << std::endl;
#pragma omp parallel for schedule(dynamic, 10000) 
  for (int i = 0; i < missingEdges_.size(); ++i) {
    ScanId row = missingEdges_[i].row;
    ScanId col = missingEdges_[i].col;
    ScanId r,c;
    if (row.fileIdx >= mergeOffset_) {
      r = mergeRoots[row];
    } else {
      r = row;
    }
    boost::unordered_map<ScanId, ScanId>::iterator it;
    if ((it = clusterMemberships.find(r)) != clusterMemberships.end()) {
      r = it->second;
    }
    if (col.fileIdx >= mergeOffset_) {
      c = mergeRoots[col];
    } else {
      c = col;
    }
    if ((it = clusterMemberships.find(c)) != clusterMemberships.end()) {
      c = it->second;
    }
    if (!(r == row && c == col)) {
      ScanId s1 = (std::min)(r,c);
      ScanId s2 = (std::max)(r,c);
      
      missingEdges_[i].row = s1;
      missingEdges_[i].col = s2;
    }
  }
}

void SparseClustering::addNewEdges(
    boost::unordered_map<ScanId, ScanId>& clusterMemberships) {
  std::cerr << "  Loading new edges." << std::endl;
  
  std::vector<PvalueTriplet> pvec;
  pvec.reserve(edgeLoadingBatchSize_);
  matrixLoader_.nextNEdges(edgeLoadingBatchSize_, pvec);
  
  size_t numNewEdges = pvec.size();
  size_t insertOffset = missingEdges_.size();
  numTotalEdges_ += numNewEdges;
#ifndef SINGLE_LINKAGE
  missingEdges_.resize(missingEdges_.size() + numNewEdges);
#endif
#pragma omp parallel for schedule(dynamic, 10000) 
  for (unsigned int i = 0; i < numNewEdges; ++i) {
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
    if (!(row == col)) {
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
  std::vector<unsigned int> limits;
  limits.push_back(0);
  std::cerr << "  Sorting " << missingEdges_.size() << " edges." << std::endl;
  std::sort(missingEdges_.begin(), missingEdges_.end(), lowerEdge);
  // make sure that identical edges end up in the same limit bin
  for (unsigned int i = 100000; i < missingEdges_.size(); i += 100000) {
    while (missingEdges_[i].col == missingEdges_[i+1].col && missingEdges_[i].row == missingEdges_[i+1].row) {
      ++i;
    }
    if (i < missingEdges_.size()) limits.push_back(i+1);
  }
  limits.push_back(missingEdges_.size());
  
  std::cerr << "  Adding new edges." << std::endl;
  // indices indicating where the last collapsed edge was inserted in missingEdges_ for each limit bin
  std::vector<size_t> collapsingIndices; 
#pragma omp parallel for schedule(dynamic, 1) 
  for (unsigned int i = 0; i < limits.size() - 1; ++i) {    
    SparseMissingEdge lastEdge = missingEdges_[limits[i]];
    lastEdge.numEdges = 0;
    size_t tmpIdx = limits[i];
    for (unsigned j = limits[i]; j < limits[i+1]; ++j) {
      if (lastEdge.col == missingEdges_[j].col && lastEdge.row == missingEdges_[j].row) {
        lastEdge.numEdges += missingEdges_[j].numEdges;
      } else {
        addNewEdge(lastEdge, tmpIdx);
        lastEdge = missingEdges_[j];
      }
    }
    addNewEdge(lastEdge, tmpIdx);
  #pragma omp critical (tmp_edges_idx_insert)
    {
      collapsingIndices.push_back(tmpIdx);
    }
  }
  
  std::cerr << "  Collapsing missing edges." << std::endl;
  std::sort(collapsingIndices.begin(), collapsingIndices.end());
  size_t curIdx = collapsingIndices.at(0);
  for (size_t i = 1; i < collapsingIndices.size(); ++i) {
    if (limits.at(i) != curIdx && limits.at(i+1) != collapsingIndices.at(i)) {
    // check if we will not accidentally overwrite data that should be moved first, otherwise execute sequentially
    #pragma omp parallel for schedule(dynamic, 1000) if (limits.at(i) - curIdx > collapsingIndices.at(i) - limits.at(i))
      for (size_t j = limits.at(i); j < collapsingIndices.at(i); ++j) {
        missingEdges_[curIdx + j - limits.at(i)] = missingEdges_.at(j);
      }
    }
    curIdx += collapsingIndices.at(i) - limits.at(i);
  }
  missingEdges_.resize(curIdx);
}

void SparseClustering::addNewEdge(const SparseMissingEdge& edge, size_t& idx) {
  if (!(edge.row == edge.col)) {
    if (edge.numEdges == clusters_[edge.row].size()*clusters_[edge.col].size()) {
      SparseEdge se(edge.value, edge.row, edge.col);
    #pragma omp critical (pq_edges_insert)
      {
        matrix_[edge.row][edge.col] = edge.value;
        matrix_[edge.col][edge.row] = edge.value;
        edgeList_.push(se);
      }
    }
#ifndef SINGLE_LINKAGE
    missingEdges_[idx++] = edge;
#endif
  }
}

void SparseClustering::updateEntry(const ScanId row, const ScanId col, 
                                   const double value) {
  edgeList_.push(SparseEdge(value, row, col));
  matrix_[row][col] = value;
  matrix_[col][row] = value;
}

// Based on http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2718652/
void SparseClustering::doClustering(double cutoff) {  
  std::cerr << "Starting MinHeap clustering" << std::endl;
  
  unsigned int itNr = 0;
  unsigned int mergeCnt = 0;
  
  typedef std::pair<ScanId, double> ColValPair;
  boost::unordered_map<ScanId, ScanId> mergeRoots;
  
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
  
  if (matrixLoader_.hasEdgesAvailable()) {
    loadNextEdges(mergeRoots);
  } else {
    std::cerr << "Could not read edges from input file." << std::endl;
    return;
  }
  while (edgeList_.top().value < cutoff) {
    SparseEdge minEdge = edgeList_.top();
    while (!isAlive_[minEdge.row] || !isAlive_[minEdge.col]) {
      if (edgeList_.size() > 0) {
        edgeList_.pop();
        minEdge = edgeList_.top();
        if (minEdge.value >= cutoff) break;
      } else {
        break;
      }
    }
    if (minEdge.value >= cutoff) break;
    if (edgeList_.size() == 0) {
      if (matrixLoader_.hasEdgesAvailable()) {
        loadNextEdges(mergeRoots);
        continue;
      } else {
        break;
      }
    }
    
    ScanId minRow = minEdge.row;
    ScanId minCol = minEdge.col;
    
    if (itNr % 10000 == 0) {
      std::cerr << "It. " << itNr << ": minRow = " << minRow << ", minCol = " << minCol << ", minEl = " << minEdge.value << std::endl;
    }
    
    ScanId minColRoot = minCol;
    if (minColRoot.fileIdx >= mergeOffset_) {
      minColRoot = mergeRoots[minCol];
    }
    
    ScanId minRowRoot = minRow;
    if (minRowRoot.fileIdx >= mergeOffset_) {
      minRowRoot = mergeRoots[minRow];
    }
    ScanId mergeScanId(mergeOffset_, mergeCnt);
    mergeRoots[mergeScanId] = minRowRoot;
    
    if (writeTree) {
      PvalueTriplet tmp(minRowRoot, minColRoot, minEdge.value);
      resultFNStream << tmp << "\n";
    }
    
    // add scans from second cluster to first one
    clusters_[mergeScanId].swap(clusters_[minRow]);
    clusters_[mergeScanId].insert(clusters_[mergeScanId].end(), 
        clusters_[minCol].begin(), clusters_[minCol].end());
    
    BOOST_FOREACH (const ColValPair & colValPair, matrix_[minRow]) {
      ScanId col = colValPair.first;
      if (col == minCol || col == minRow || !isAlive_[col]) continue;
      if (matrix_[minCol].find(col) != matrix_[minCol].end()) {
        //double value = (colValPair.second * n + matrix_[minCol][col] * m)/(n+m); // UPGMA
#ifdef SINGLE_LINKAGE
        double value = (std::min)(colValPair.second, matrix_[minCol][col]); // single linkage
#else
        double value = (std::max)(colValPair.second, matrix_[minCol][col]); // complete linkage
#endif
        updateEntry(mergeScanId, col, value);
      }
    }
    
    edgeList_.pop();
    isAlive_[minRow] = false;
    isAlive_[minCol] = false;
    isAlive_[mergeScanId] = true;
    
    // TODO: make this a batch clean up (say every 1000 iterations) to increase speed
    clusters_.erase(minRow);
    clusters_.erase(minCol);
    matrix_.erase(minRow);
    matrix_.erase(minCol);
    
    ++mergeCnt;
    ++itNr;
  }
  
  std::cerr << "Finished MinHeap clustering" << std::endl;
  
  if (writeMissingEdges_) {
    loadNextEdges(mergeRoots);
    std::sort(missingEdges_.begin(), missingEdges_.end(), lowerPval);
    std::string clusterPairMissingFN = clusterPairFN_ + ".missing_edges.tsv";
    std::ofstream resultMissingFNStream(clusterPairMissingFN.c_str());
    if (resultMissingFNStream.is_open()) {
      BOOST_FOREACH (SparseMissingEdge& sme, missingEdges_) {
        if (sme.value < cutoff - 5.0) {
          if (sme.row.fileIdx >= mergeOffset_) {
            sme.row = mergeRoots[sme.row];
          }
          if (sme.col.fileIdx >= mergeOffset_) {
            sme.col = mergeRoots[sme.col];
          }   
          resultMissingFNStream << sme.row << "\t" << sme.col << "\t" << sme.value << "\t" << sme.numEdges << "\n";
        } else {
          break;
        }
      }
    }
  }
  
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

void SparseClustering::readMatrix(std::string & matrixInFN) {  
  std::ifstream inputStream(matrixInFN.c_str());
  std::string line, colValString;
  while (std::getline(inputStream, line)) {
    double logPval;
    bool first = true;
    std::istringstream ss(line);
    unsigned int row;
    ss >> row;
    while (ss >> colValString) { // every value after a colon is a p-value followed by a space
      std::size_t found = colValString.find_first_of(":");
      unsigned int col = boost::lexical_cast<unsigned int>(colValString.substr(0, found));
      double logPval = boost::lexical_cast<double>(colValString.substr(found+1));
      if (row != col) {
        // TODO: Fix this!
        minInsert(ScanId(0,row), ScanId(0, col), logPval);
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
  
  matrix.minInsert(ScanId(0,1),ScanId(0,1),-20.0);
  matrix.minInsert(ScanId(0,1),ScanId(0,2),-2.0);
  matrix.minInsert(ScanId(0,1),ScanId(0,3),-6.0);
  matrix.minInsert(ScanId(0,1),ScanId(0,4),-6.0);
  matrix.minInsert(ScanId(0,2),ScanId(0,3),-1.0);
  matrix.minInsert(ScanId(0,2),ScanId(0,4),-10.0);
  matrix.minInsert(ScanId(0,3),ScanId(0,4),-4.0);
  
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
  if (clusterRow[0] == ScanId(0,1) 
      && clusterRow[1] == ScanId(0,3) 
      && clusterRow[2] == ScanId(0,2) 
      && clusterRow[3] == ScanId(0,4)) {
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


