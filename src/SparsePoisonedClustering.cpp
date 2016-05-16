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
 
#include "SparsePoisonedClustering.h"

void SparsePoisonedClustering::loadNextEdges() {
  edgesLeft_ = false;
  std::cerr << "  Loading new edges." << std::endl;
  
  size_t numNewEdges = pvals_.size();
  numTotalEdges_ += numNewEdges;
  for (int i = 0; i < numNewEdges; ++i) {
    ScanId row = pvals_[i].scannr1;
    ScanId col = pvals_[i].scannr2;
    
    isAlive_[row] = true;
    isAlive_[col] = true;
    
    ScanId s1 = (std::min)(row, col);
    ScanId s2 = (std::max)(row, col);
#ifdef SINGLE_LINKAGE
    if (row != col) {
      size_t idx = 0;
      addNewEdge(SparseMissingEdge(pvec[i].pval, s1, s2, clusters_[s1].size()*clusters_[s2].size()), idx);
    }
#else
    insertEdge(s1, s2, pvals_[i].pval);
#endif
  }
  matrix_.sortRows();
  pvals_.clear();
}

void SparsePoisonedClustering::doClustering(double cutoff) {  
  std::cerr << "Starting MinHeap clustering" << std::endl;
  
  // decide if we write the results to a file or stdout
  std::ofstream resultFNStream;
  bool writeTree = false;
  if (clusterPairFN_.size() > 0) {
    resultFNStream.open(clusterPairFN_.c_str(), std::ios_base::app);
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
      if (isPoisoned_[minEdge.row] || isPoisoned_[minEdge.col]) {
        isPoisoned_[minEdge.row] = true;
        isPoisoned_[minEdge.col] = true;
        
        ScanId minRowRoot = getRoot(minEdge.row);
        ScanId minColRoot = getRoot(minEdge.col);
        poisonedEdges_.push_back(PvalueTriplet(std::min(minRowRoot, minColRoot), 
                                               std::max(minRowRoot, minColRoot), 
                                               minEdge.value));
      } else {
        isAlive_[minEdge.row] = false;
        isAlive_[minEdge.col] = false;
        
        if (mergeCnt % 10000 == 0) {
          std::cerr << "It. " << mergeCnt << ": minRow = " << minEdge.row 
                    << ", minCol = " << minEdge.col 
                    << ", minEl = " << minEdge.value << std::endl;
        }
        
        ScanId minRowRoot = getRoot(minEdge.row);
        ScanId mergeScanId(mergeOffset_, mergeCnt++);
        mergeRoots_[mergeScanId] = minRowRoot;
        isAlive_[mergeScanId] = true;
        
        if (writeTree) {
          ScanId minColRoot = getRoot(minEdge.col);
          PvalueTriplet tmp(minRowRoot, minColRoot, minEdge.value);
          resultFNStream << tmp << "\n";
        }
        
        updateMatrix(minEdge.row, minEdge.col, mergeScanId);
      }
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
