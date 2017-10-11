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
 
#ifndef MARACLUSTER_RANKMERGE_H_
#define MARACLUSTER_RANKMERGE_H_

#include <iostream>
#include <string>
#include <vector>

#include <boost/foreach.hpp>

#include "MZIntensityPair.h"
#include "SpectrumHandler.h"
#include "BinAndRank.h"

namespace maracluster {

class RankMerge {
 public:
  inline static void setDpThresh(double _dpThresh) { dpThresh = _dpThresh; }
  static void merge(std::vector<MZIntensityPair>& mziPairsIn,
                    std::vector<MZIntensityPair>& mziPairsFrom, 
                    BinRanks& rankIn, BinRanks& rankFrom, double weight,
                    double shift = 0.0, double scaling = 1.0);
  static void mergeMinMax(std::vector<MZIntensityPair>& mziPairsIn,
                          std::vector<MZIntensityPair>& mziPairsFrom, 
                          double weight, bool regionRanks, 
                          double shift = 0.0, double scaling = 1.0);
 protected:
  static double dpThresh;
  static int maxRankDiff;
  static inline bool ranksInRange(unsigned int binRankIn, 
                                  unsigned int binRankFrom) {
    return (binRankIn > 0 && binRankFrom > 0 && 
            std::abs((int)binRankIn - (int)binRankFrom) <= maxRankDiff);
  }
};

} /* namespace maracluster */

#endif /* MARACLUSTER_RANKMERGE_H_ */
