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
 
#ifndef SPARSE_POISONED_CLUSTERING_H
#define SPARSE_POISONED_CLUSTERING_H

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

#include "SparseClustering.h"

class SparsePoisonedClustering : public SparseClustering {
 public:
  SparsePoisonedClustering() : SparseClustering(), edgesLeft_(false) { }
  
  void markPoisoned(const ScanId& scanId) {
    isPoisoned_[scanId] = true;
  }
  
  void initPvals(std::vector<PvalueTriplet>& pvec) {
    edgesLeft_ = true;
    pvals_.swap(pvec);
  }
  
  void getPoisonedEdges(std::vector<PvalueTriplet>& pvec) {
    pvec.swap(poisonedEdges_);
  }
  
  void doClustering(double cutoff);
  
 protected:  
  boost::unordered_map<ScanId, bool> isPoisoned_;
  std::vector<PvalueTriplet> pvals_, poisonedEdges_;
  bool edgesLeft_;
  
  void loadNextEdges();
  bool edgesLeft() { return edgesLeft_; }
};

#endif // SPARSE_POISONED_CLUSTERING_H

