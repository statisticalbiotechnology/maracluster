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
 
#ifndef MARACLUSTER_MSFILEHANDLER_H_
#define MARACLUSTER_MSFILEHANDLER_H_

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <cstdio>

#include "pwiz/data/msdata/MSDataFile.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include "SpectrumHandler.h"
#include "MZIntensityPair.h"
#include "BinAndRank.h"
#include "BinSpectra.h"
#include "PeakCounts.h"
#include "SpectrumFileList.h"
#include "ScanMergeInfo.h"
#include "PvalueCalculator.h"

#include "BatchSpectrum.h"

namespace maracluster {

class MSFileHandler {
 public:
  static unsigned int chargeErrorTolerance_;
  static bool splitMassChargeStates_;
  
  MSFileHandler(std::string& spectrumOutFN) : spectrumOutFN_(spectrumOutFN) {}
  
  /* Splits mass charge states for MSGF+ input */
  void msgfFixMzML(const std::string& spectrumInFN);
  
  void calcRankDotProducts(const std::string& spectrumInFN);
  
  inline void setCombineSets(SpectrumFileList& fileList, 
                             std::vector<ScanMergeInfoSet>& combineSets) {
    fileList_ = fileList;
    combineSets_ = combineSets;
  }
  
  static bool validMs2OutputFN(std::string& outputFN);
  static std::string getOutputFormat(const std::string& outputFN);
  static void calcPeakCount(pwiz::msdata::SpectrumListPtr specList, 
      PeakCounts& mzMap,
      std::string resultFN = "", bool normalizeXCorr = false, 
      bool addPrecursorMass = true, bool truncatePeaks = true, 
      unsigned int chargeFilter = 0);
 protected:
  SpectrumFileList fileList_;
  std::vector<ScanMergeInfoSet> combineSets_;
  std::string spectrumOutFN_;
  
  static size_t getSpectrumIdxFromScannr(pwiz::msdata::SpectrumListPtr sl, 
                                         unsigned int scannr);
  static void writeMSData(pwiz::msdata::MSData& msd, const std::string& outputFN);
    
  void addSpectrumWithMccs(pwiz::msdata::SpectrumPtr consensusSpec, 
     std::vector<MassChargeCandidate>& consensusMccs,
     unsigned int scannr, pwiz::msdata::SpectrumListSimplePtr mergedSpectra);
};

} /* namespace maracluster */

#endif /* MARACLUSTER_MSFILEHANDLER_H_ */
