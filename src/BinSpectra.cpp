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
 
#include "BinSpectra.h"

namespace maracluster {

const double BinSpectra::kBinWidth = 1.000508;
const double BinSpectra::kBinShift = 0.32;

// const double BinSpectra::kBinWidth = 0.05;
// const double BinSpectra::kBinShift = 0.0;

const unsigned int BinSpectra::kRankWindow = 10u;
const unsigned int BinSpectra::kMaxRank = 3u;

// this function requires that mziPairs are sorted by mz value
unsigned int BinSpectra::bin(std::vector<MZIntensityPair>& mziPairsIn, std::vector<BinnedMZIntensityPair>& mziPairsBinned) {
  unsigned int lastBin = 0;
  double acc = 0.0;
  BOOST_FOREACH (const MZIntensityPair& mziPair, mziPairsIn) {
    unsigned int currentBin = getBin(mziPair.mz);
    if (currentBin == lastBin) {
      acc = (std::max)(mziPair.intensity, acc);
    } else {
      if (lastBin > 0) {
        mziPairsBinned.push_back(BinnedMZIntensityPair(lastBin,acc));
      }
      acc = mziPair.intensity;
      lastBin = currentBin;
    }
  }
  return lastBin;
}

unsigned int BinSpectra::binDense(std::vector<MZIntensityPair>& mziPairsIn, std::vector<BinnedMZIntensityPair>& mziPairsBinned) {
  unsigned int lastBin = 0;
  double acc = 0.0;
  BOOST_FOREACH (const MZIntensityPair& mziPair, mziPairsIn) {
    unsigned int currentBin = getBin(mziPair.mz);
    
    unsigned int emptyBin = lastBin + 1;
    while (emptyBin < currentBin) {
       mziPairsBinned.push_back(BinnedMZIntensityPair(emptyBin,0));
       ++emptyBin;
    }
    if (currentBin == lastBin) {
      acc = (std::max)(mziPair.intensity, acc);
    } else {
      if (lastBin > 0) {
        mziPairsBinned.push_back(BinnedMZIntensityPair(lastBin,acc));
      }
      acc = mziPair.intensity;
      lastBin = currentBin;
    }
  }
  return lastBin;
}

void BinSpectra::binBinary(std::vector<MZIntensityPair>& mziPairsIn, std::vector<unsigned int>& peakBins, double intThresh) {
  unsigned int lastBin = 0;
  BOOST_FOREACH (const MZIntensityPair& mziPair, mziPairsIn) {
    if (mziPair.intensity > intThresh) {
      unsigned int currentBin = getBin(mziPair.mz);
      if (currentBin != lastBin) {
        peakBins.push_back(currentBin);
        lastBin = currentBin;
      }
    }
  }
}

#ifdef DOT_PRODUCT
void BinSpectra::binBinaryTruncated(std::vector<MZIntensityPair>& mziPairs, std::vector<unsigned int>& peakBins, 
                                      const unsigned int nPeaks, double precMass) {
  std::map<unsigned int, bool> peakFound;
  peakBins.clear();
  std::sort( mziPairs.begin(), mziPairs.end(), SpectrumHandler::greaterIntensity );
  unsigned int peakCnt = 0;
  double maxIntensity = mziPairs.begin()->intensity;
  BOOST_FOREACH (const MZIntensityPair& mziPair, mziPairs) {
    if (mziPair.mz < precMass) {
      unsigned int bin = getBin(mziPair.mz);
      if (!peakFound[bin]) {
        peakFound[bin] = true;
        peakBins.push_back(bin);
        peakBins.push_back(static_cast<unsigned int>(std::sqrt(mziPair.intensity/maxIntensity)*1e4));
        ++peakCnt;
      }
      if (peakCnt >= nPeaks) break;
    }
  }
}
#else
void BinSpectra::binBinaryTruncated(std::vector<MZIntensityPair>& mziPairs, std::vector<unsigned int>& peakBins, 
                                      const unsigned int nPeaks, double precMass) {
  std::map<unsigned int, bool> peakFound;
  peakBins.clear();
  std::sort( mziPairs.begin(), mziPairs.end(), SpectrumHandler::greaterIntensity );
  unsigned int peakCnt = 0;
  BOOST_FOREACH (const MZIntensityPair& mziPair, mziPairs) {
    if (mziPair.mz < precMass) {
      unsigned int bin = getBin(mziPair.mz);
      if (!peakFound[bin]) {
        peakFound[bin] = true;
        peakBins.push_back(bin);
        ++peakCnt;
      }
      if (peakCnt >= nPeaks) break;
    }
  }
  std::sort(peakBins.begin(), peakBins.end());
}
#endif

void BinSpectra::binBinaryPeakPicked(std::vector<MZIntensityPair>& mziPairs, std::vector<unsigned int>& peakBins, 
                                      const unsigned int nPeaks, double precMass, bool reportDuplicates) {
  std::map<unsigned int, bool> peakFound;
  std::map<unsigned int, unsigned int> peakRank;
  std::sort( mziPairs.begin(), mziPairs.end(), SpectrumHandler::greaterIntensity );
  peakBins.clear();
  unsigned int peakCnt = 0;
  BOOST_FOREACH (const MZIntensityPair& mziPair, mziPairs) {
    if (mziPair.mz < precMass) {
      unsigned int bin = getBin(mziPair.mz);
      
      unsigned int startWindowBin = 0;
      if (bin > kRankWindow/2) startWindowBin = bin - kRankWindow/2;
      unsigned int endWindowBin = bin + kRankWindow/2;
      for (unsigned int windowBin = startWindowBin; windowBin <= endWindowBin; ++windowBin) {
        peakRank[windowBin] += 1;
      }
      
      if (peakRank[bin] <= kMaxRank) {
        if (!peakFound[bin] || reportDuplicates) {
          peakFound[bin] = true;
          peakBins.push_back(bin);
        }
        ++peakCnt;
        if (peakCnt >= nPeaks) break; 
      }
    }
  }
  std::sort(peakBins.begin(), peakBins.end());
}

void BinSpectra::printIntensities(std::vector<BinnedMZIntensityPair>& mziPairsBinned) {
  std::vector<double> intensities(2000);
  BOOST_FOREACH (BinnedMZIntensityPair& mziPair, mziPairsBinned) {
		intensities[mziPair.binIdx] = mziPair.intensity;
	}
	BOOST_FOREACH (double intensity, intensities) {
	  std::cout << intensity << "\t";
	}
	std::cout << std::endl;
}

} /* namespace maracluster */
