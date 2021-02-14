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

namespace maracluster {

bool SparsePoisonedClustering::verifyEdgeBoth(SparseEdge& minEdge) {
  if (isPoisoned(minEdge.row) && isPoisoned(minEdge.col)) {
    ScanId minRowRoot = getRoot(minEdge.row);
    ScanId minColRoot = getRoot(minEdge.col);
    poisonedEdges_.push_back(PvalueTriplet(std::min(minRowRoot, minColRoot), 
                                           std::max(minRowRoot, minColRoot), 
                                           minEdge.value));
    
    if (poisonedEdges_.size() % 10000000 == 0) {
      std::cerr << "Collected " << poisonedEdges_.size() << " poisoned edges" << std::endl;
    }
    
    return false;
  }
  return true;
}

bool SparsePoisonedClustering::verifyEdgeSingle(SparseEdge& minEdge) {
  if (isPoisoned(minEdge.row) || isPoisoned(minEdge.col)) {
    markPoisoned(minEdge.row);
    markPoisoned(minEdge.col);
    
    ScanId minRowRoot = getRoot(minEdge.row);
    ScanId minColRoot = getRoot(minEdge.col);
    poisonedEdges_.push_back(PvalueTriplet(std::min(minRowRoot, minColRoot), 
                                           std::max(minRowRoot, minColRoot), 
                                           minEdge.value));
    
    if (poisonedEdges_.size() % 10000000 == 0) {
      std::cerr << "Collected " << poisonedEdges_.size() << " poisoned edges" << std::endl;
    }
    
    return false;
  }
  return true;
}

} /* namespace maracluster */
