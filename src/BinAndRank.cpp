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
 
#include "BinAndRank.h"

namespace maracluster {

const unsigned int BinAndRank::kNumRegions = 10u;

void BinAndRank::binAndXCorrNormalize(std::vector<MZIntensityPair>& mziPairsIn, 
                                        std::vector<BinnedMZIntensityPair>& mziPairsNormalized) {
  std::vector<BinnedMZIntensityPair> mziPairsBinned;
  BinSpectra::binDense(mziPairsIn, mziPairsBinned);
  
  BOOST_FOREACH (MZIntensityPair& mziPair, mziPairsBinned) {
    mziPair.intensity = sqrt(mziPair.intensity);
  }
  
  // find maximum in each of the 10 regions
  unsigned int numRegions = 10;
  unsigned int numBins = mziPairsBinned.size();
  unsigned int regionSize = numBins/numRegions;
  std::vector<double> regionMax(numRegions);
  for (unsigned int region = 0; region < numRegions; ++region) {
    regionMax[region] = std::max_element(mziPairsBinned.begin() + region*regionSize, mziPairsBinned.begin() + (region+1)*regionSize, SpectrumHandler::lessIntensity)->intensity;
    if (regionMax[region] == 0.0) regionMax[region] = 1.0;
  }
  
  // normalize to regionMax and discard peaks below 0.05*globalMax
  double globalMax = *std::max_element(regionMax.begin(), regionMax.end());
  BOOST_FOREACH (BinnedMZIntensityPair& mziPair, mziPairsBinned) {
    if (mziPair.intensity / globalMax < 0.05) {
      mziPair.intensity = 0;
    } else {
      mziPair.intensity /= regionMax[getRegion(mziPair.binIdx, mziPairsBinned.size())];
    }
  }
  
  // apply smoothing with a 151 bin window
  mziPairsNormalized.clear();
  double runningAvg = 0.0;
  for (std::vector<BinnedMZIntensityPair>::iterator it = mziPairsBinned.begin(); it != mziPairsBinned.begin() + 75; ++it) {
    runningAvg += it->intensity/151;
  }
  for (std::vector<BinnedMZIntensityPair>::iterator it = mziPairsBinned.begin(); it != mziPairsBinned.end(); ++it) {
    mziPairsNormalized.push_back(BinnedMZIntensityPair(it->binIdx, it->intensity - runningAvg));
    if (it + 75 < mziPairsBinned.end()) {
      runningAvg += (it+75)->intensity/151;
    }
    if (it > mziPairsBinned.begin() + 75) {
      runningAvg -= (it-75)->intensity/151;
    }
  }
  
  // keep only the 100 bins with the highest intensity
  std::sort( mziPairsNormalized.begin(), mziPairsNormalized.end(), SpectrumHandler::greaterIntensity );
  mziPairsNormalized.resize((std::min)(100u, static_cast<unsigned int>(mziPairsNormalized.size())));
}


void BinAndRank::binAndRank(std::vector<MZIntensityPair>& mziPairsIn, BinRanks& mzRankMap) {
  std::vector<BinnedMZIntensityPair> mziPairsBinned;
  BinSpectra::bin(mziPairsIn, mziPairsBinned);
  std::sort( mziPairsBinned.begin(), mziPairsBinned.end(), SpectrumHandler::greaterIntensity );
  unsigned int rank = 1;
  BOOST_FOREACH (const BinnedMZIntensityPair& mziPair, mziPairsBinned) {
    mzRankMap[mziPair.binIdx] = rank;
    ++rank;
  }
}

void BinAndRank::binAndRegionRank(std::vector<MZIntensityPair>& mziPairsIn, BinRanks& mzRankMap) {
  std::vector<BinnedMZIntensityPair> mziPairsBinned;
  unsigned int maxBin = BinSpectra::bin(mziPairsIn, mziPairsBinned);
  std::sort( mziPairsBinned.begin(), mziPairsBinned.end(), SpectrumHandler::greaterIntensity );
  std::vector<unsigned int> rank(10);
  BOOST_FOREACH (const BinnedMZIntensityPair& mziPair, mziPairsBinned) {
  //TODO: fix maxBin to only the first spectrum, otherwise the ranks spill over the regions
    mzRankMap[mziPair.binIdx] = ++rank[getRegion(mziPair.mz, maxBin)];
  }
}

void BinAndRank::binAndRank(std::vector<MZIntensityPair>& mziPairsIn, BinRankScores& mzRankMap) {
  std::vector<BinnedMZIntensityPair> mziPairsBinned;
  BinSpectra::bin(mziPairsIn, mziPairsBinned);
  std::sort( mziPairsBinned.begin(), mziPairsBinned.end(), SpectrumHandler::greaterIntensity );
  unsigned int rank = 1;
  BOOST_FOREACH (const BinnedMZIntensityPair& mziPair, mziPairsBinned) {
    mzRankMap[mziPair.binIdx] = rankToScore(rank);
    ++rank;
    if (rank > 50) break;
  }
}

void BinAndRank::normalizeRankScores(BinRankScores& rankScores) {
  BOOST_FOREACH(BinRankScores::value_type& mzRankPair, rankScores) {
    for (int i = (std::max)(0,(int)mzRankPair.first - 75); i <= (std::min)(2000,(int)mzRankPair.first + 75); ++i) {
      mzRankPair.second -= rankScores[i]/151;
    }
  }
}

double BinAndRank::rankDotProduct(BinRanks& rankIn, BinRanks& rankFrom) {
  double dotProduct = 0;
  BOOST_FOREACH(BinRanks::value_type& mzRankPair, rankIn) {
    dotProduct += rankToScore(mzRankPair.second) * rankToScore(rankFrom[mzRankPair.first]);
  }
  return dotProduct;
}

double BinAndRank::rankDotProduct(BinRankScores& rankIn, BinRankScores& rankFrom) {
  double dotProduct = 0;
  BOOST_FOREACH(BinRankScores::value_type& mzRankPair, rankIn) {
    dotProduct += mzRankPair.second * rankFrom[mzRankPair.first];
  }
  return dotProduct;
}

std::ostream &operator<<(std::ostream &out, std::vector<MZIntensityPair>& mziPairs) {
  BOOST_FOREACH(const MZIntensityPair & mziPair, mziPairs) {
    out << mziPair.mz << '\t' << mziPair.intensity << std::endl;
  }
  return out;
}

std::ostream &operator<<(std::ostream &out, BinRanks &ranks) {
  for (unsigned int i = 0; i < 2000; ++i) {
    out << ranks[i] << '\t';
  }
  return out;
}

std::ostream &operator<<(std::ostream &out, BinRankScores &rankScores) {
  for (unsigned int i = 0; i < 2000; ++i) {
    out << rankScores[i] << '\t';
  }
  return out;
}

} /* namespace maracluster */
