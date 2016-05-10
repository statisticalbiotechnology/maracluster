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
 
#ifndef BIN_AND_RANK_H
#define BIN_AND_RANK_H

#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ostream>

#include <boost/foreach.hpp>

#include "SpectrumHandler.h"
#include "BinSpectra.h"
#include "MZIntensityPair.h"

typedef std::map<unsigned int, unsigned int> BinRanks;
typedef std::map<unsigned int, double> BinRankScores;

class BinAndRank {

  public:    
    static void binAndXCorrNormalize(std::vector<MZIntensityPair>& mziPairsIn, std::vector<BinnedMZIntensityPair>& mziPairsNormalized);
    
    static void binAndRank(std::vector<MZIntensityPair>& mziPairsIn, BinRanks& mzRankMap);
    static void binAndRegionRank(std::vector<MZIntensityPair>& mziPairsIn, BinRanks& mzRankMap);
    
    static void binAndRank(std::vector<MZIntensityPair>& mziPairsIn, BinRankScores& mzRankMap);
    static void normalizeRankScores(BinRankScores& rankScores);
    
    static double rankDotProduct(BinRanks& rankIn, BinRanks& rankFrom);
    static double rankDotProduct(BinRankScores& rankIn, BinRankScores& rankFrom);
    
  protected:
    static const unsigned int kNumRegions;
    static inline unsigned int getRegion(double mz, unsigned int maxBin) {
      return (std::min)(kNumRegions-1, BinSpectra::getBin(mz)/(maxBin/kNumRegions));
    }
    static inline unsigned int getBinRegion(unsigned int bin, unsigned int maxBin) { 
      return (std::min)(kNumRegions-1, bin/(maxBin/kNumRegions)); 
    }
    static inline double rankToScore(unsigned int rank) { 
      if (rank > 0) 
        return (std::max)(0.0,(51.0 - rank)/50.0);
      else
        return 0;
    }
};

std::ostream &operator<<(std::ostream &out, std::vector<MZIntensityPair>& mziPairs);
std::ostream &operator<<(std::ostream &out, BinRanks &ranks);
std::ostream &operator<<(std::ostream &out, BinRankScores &rankScores);

#endif // BIN_AND_RANK_H
