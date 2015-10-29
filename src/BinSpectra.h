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
 
#ifndef BIN_SPECTRA_H
#define BIN_SPECTRA_H

#include <vector>
#include <iostream>
#include <cstring>
#include <boost/foreach.hpp>
#include "MZIntensityPair.h"
#include "SpectrumHandler.h"

struct BinnedMZIntensityPair; // forward declaration

class BinSpectra {
  public:
    static const double kBinWidth, kBinShift;
    
    BinSpectra() {}
    
    static inline unsigned int getBin(double mz, double binWidth = kBinWidth, double binShift = kBinShift) { 
      return static_cast<unsigned int>((mz / binWidth) + binShift);
    }

    static inline double getMZ(unsigned int bin, double binWidth = kBinWidth, double binShift = kBinShift) {
      return (bin + 0.5 - binShift)*binWidth;
    }
    
    static unsigned int bin(std::vector<MZIntensityPair>& mziPairsIn, std::vector<BinnedMZIntensityPair>& mziPairsBinned);
    static unsigned int binDense(std::vector<MZIntensityPair>& mziPairsIn, std::vector<BinnedMZIntensityPair>& mziPairsBinned);
    static void binBinary(std::vector<MZIntensityPair>& mziPairsIn, std::vector<unsigned int>& peakBins, double intThresh = 0.0);
    static void binBinaryTruncated(std::vector<MZIntensityPair>& mziPairsIn, std::vector<unsigned int>& peakBins, 
                                      const unsigned int nPeaks, double precMass);
    static void binBinaryPeakPicked(std::vector<MZIntensityPair>& mziPairsIn, std::vector<unsigned int>& peakBins, 
                                      const unsigned int nPeaks, double precMass, bool reportDuplicates = false);
    static void printIntensities(std::vector<BinnedMZIntensityPair>& mziPairsBinned);
  protected:
    static const unsigned int kRankWindow, kMaxRank;
};

struct BinnedMZIntensityPair : public MZIntensityPair {
  public:
    unsigned int binIdx;
  
    BinnedMZIntensityPair(unsigned int binIdx_, double mz_, double intensity_, double multiplicity_) {
      init(binIdx_, mz_, intensity_, multiplicity_);
    }
    
    BinnedMZIntensityPair(unsigned int binIdx_, double mz_, double intensity_) {
      init(binIdx_, mz_, intensity_, 1.0);
    }
    
    BinnedMZIntensityPair(unsigned int binIdx_, double intensity_) {
      init(binIdx_, BinSpectra::getMZ(binIdx_), intensity_, 1.0);
    }
    
    BinnedMZIntensityPair() {
      init(0,0.0,0.0,1.0);
    }
    
  protected:
    void init(unsigned int binIdx_, double mz_, double intensity_, double multiplicity_) {
      binIdx = binIdx_;
      mz = mz_;
      intensity = intensity_;
      multiplicity = multiplicity_;
    }
};


#endif // BIN_SPECTRA_H
