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

#ifndef MARACLUSTER_SPARSEPOISONEDCLUSTERING_H_
#define MARACLUSTER_SPARSEPOISONEDCLUSTERING_H_

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

#include "SparseClustering.h"

namespace maracluster {

/**
  This class extends the SparseClustering class by adding functionality for 
  "poisoned" clustering. A set of nodes (ScanIds) are marked as poisoned, 
  e.g. since their precursor m/z is close to the bin/batch m/z border. Edges
  that have at least one node that is poisoned will not be clustered, and
  the other node for that edge will also be marked as poisoned (hence the
  poison analogy). These poisoned nodes and edges can then be clustered at a 
  later stage.
*/
class SparsePoisonedClustering : public SparseClustering {
 public:
  SparsePoisonedClustering() : SparseClustering() { }
  
  inline void markPoisoned(const ScanId& scanId) {
    isPoisoned_[scanId] = true;
  }
  
  inline bool isPoisoned(const ScanId& scanId) {
    return isPoisoned_[scanId];
  }
  
  void initPvals(std::vector<PvalueTriplet>& pvec) {
    matrixLoader_.initVector(pvec);
  }
  
  void getPoisonedEdges(std::vector<PvalueTriplet>& pvec) {
    pvec.swap(poisonedEdges_);
  }
  
 protected:  
  boost::unordered_map<ScanId, bool> isPoisoned_;
  std::vector<PvalueTriplet> poisonedEdges_;
  
  bool verifyEdgeBoth(SparseEdge& minEdge);
  bool verifyEdgeSingle(SparseEdge& minEdge);
};

} /* namespace maracluster */

#endif /* MARACLUSTER_SPARSEPOISONEDCLUSTERING_H_ */

