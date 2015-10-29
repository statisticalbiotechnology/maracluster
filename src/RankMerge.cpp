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
 
#include "RankMerge.h"

double RankMerge::dpThresh = 12.0;
unsigned int RankMerge::maxRankDiff = 10u;

void RankMerge::merge(std::vector<MZIntensityPair>& mziPairsIn, std::vector<MZIntensityPair>& mziPairsFrom, 
                              BinRanks& rankIn, BinRanks& rankFrom, double weight, double shift, double scaling) {
  unsigned int fromCount = 0;
  BOOST_FOREACH (MZIntensityPair& mziPair, mziPairsIn) {
    while (mziPairsFrom.at(fromCount).mz < (mziPair.mz + shift)*scaling) {
      ++fromCount;
      if (fromCount >= mziPairsFrom.size()) break;
    }
    if (fromCount >= mziPairsFrom.size()) break;
    if (fromCount >= 1) {
      unsigned int mzBin = BinSpectra::getBin(mziPair.mz);
      if (ranksInRange(rankIn[mzBin], rankFrom[mzBin])) {
        mziPair.intensity += weight * SpectrumHandler::interpolateIntensity(mziPairsFrom.at(fromCount-1), mziPairsFrom.at(fromCount), (mziPair.mz + shift)*scaling);
      }
    }
  }
}

void RankMerge::mergeMinMax(std::vector<MZIntensityPair>& mziPairsIn, std::vector<MZIntensityPair>& mziPairsFrom, 
                                    double weight, bool regionRanks, double shift, double scaling) {
  BinRanks rankIn, rankFrom;
  double dotProduct;
  if (regionRanks) {
    BinRanks dpRankIn, dpRankFrom;
    BinAndRank::binAndRank(mziPairsIn, dpRankIn);
    BinAndRank::binAndRank(mziPairsFrom, dpRankFrom);
    dotProduct = BinAndRank::rankDotProduct(dpRankIn, dpRankFrom);
    
    BinAndRank::binAndRegionRank(mziPairsIn, rankIn);
    BinAndRank::binAndRegionRank(mziPairsFrom, rankFrom);
  } else {
    BinAndRank::binAndRank(mziPairsIn, rankIn);
    BinAndRank::binAndRank(mziPairsFrom, rankFrom);
  }
  
  if (dotProduct >= dpThresh) {
    unsigned int fromCount = 0;
    BOOST_FOREACH (MZIntensityPair& mziPair, mziPairsIn) {
      while (mziPairsFrom.at(fromCount).mz < (mziPair.mz + shift)*scaling) {
        ++fromCount;
        if (fromCount >= mziPairsFrom.size()) break;
      }
      if (fromCount >= mziPairsFrom.size()) break;
      if (fromCount >= 1) {
        unsigned int mzBin = BinSpectra::getBin(mziPair.mz);
        double interpolatedIntensity = SpectrumHandler::interpolateIntensity(mziPairsFrom.at(fromCount-1), mziPairsFrom.at(fromCount), (mziPair.mz + shift)*scaling);
        if (ranksInRange(rankIn[mzBin-1], rankFrom[mzBin-1]) || ranksInRange(rankIn[mzBin], rankFrom[mzBin]) || ranksInRange(rankIn[mzBin+1], rankFrom[mzBin+1])) { 
        //if (ranksInRange(rankIn[mzBin], rankFrom[mzBin])) { 
          mziPair.intensity = (std::max)(interpolatedIntensity, mziPair.intensity);
        } else {
          mziPair.intensity = (std::min)(interpolatedIntensity, mziPair.intensity);
        }
      }
    }
  }
}
