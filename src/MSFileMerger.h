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

#ifndef MARACLUSTER_MSFILEMERGER_H_
#define MARACLUSTER_MSFILEMERGER_H_

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <cstdio>

#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/MSDataMerger.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include "Globals.h"
#include "MSFileHandler.h"
#include "MSClusterMerge.h"
#include "ScanId.h"

namespace maracluster {

struct MergeScanIndex {
  MergeScanIndex(unsigned int _spectrumIndex, unsigned int _posInCluster, unsigned int _mergeIdx) :
    spectrumIndex(_spectrumIndex), posInCluster(_posInCluster), mergeIdx(_mergeIdx) {}
  unsigned int spectrumIndex, posInCluster, mergeIdx;
};

class MSFileMerger : public MSFileHandler {
 public:
  static int mergeMethod_;
  static bool normalize_;
  static int maxMSFilePtrs_, maxSpectraPerFile_;
  static unsigned int maxConsensusSpectraPerFile_;

  MSFileMerger(std::string& spectrumOutFN) :
    numClusterBins_(0), numBatches_(0),
    numMSFilePtrsPerBatch_(0), MSFileHandler(spectrumOutFN) {}

  void parseClusterFileForMerge(const std::string& clusterFile,
      const size_t minClusterSize);

  void mergeSpectra();

  void mergeAllSpectra(const std::string& spectrumInFN);

  void mergeSpectraSmall();
  void mergeSpectraScalable();

  inline size_t getClusterBin(const ScanId& si) {
    return si.scannr % numClusterBins_;
  }

  inline static bool lessIndex(const MergeScanIndex& a,
    const MergeScanIndex& b) { return (a.spectrumIndex < b.spectrumIndex); }
 protected:
  unsigned int numClusterBins_, numBatches_, numMSFilePtrsPerBatch_;

  void mergeTwoSpectra(
      std::vector<MZIntensityPair>& mziPairsIn,
      std::vector<MZIntensityPair>& mziPairsFrom,
      double weight);

  void mergeSpectraSet(std::vector<pwiz::msdata::SpectrumPtr>& spectra,
                              pwiz::msdata::SpectrumPtr& consensusSpec);
  void mergeSpectraSetMSCluster(
      std::vector<pwiz::msdata::SpectrumPtr>& spectra,
      ScanId scannr, pwiz::msdata::SpectrumListSimplePtr mergedSpectra);

  void mergeSpectraBin(size_t clusterBin, 
    std::vector<pwiz::msdata::MSDataPtr>& msdVector,
    pwiz::msdata::SpectrumListSimplePtr mergedSpectra);

  void mergeSplitSpecFiles();
  
  inline size_t calcNumBatches(size_t total, size_t batchSize) {
    return (total - 1u) / batchSize + 1u;
  }
  virtual void mergeMccs(std::vector<MassChargeCandidate>& allMccs,
      std::vector<MassChargeCandidate>& consensusMccs, int clusterIdx);
  
  void splitSpecFilesByConsensusSpec(
    std::map<ScanId, ScanId>& scannrToMergedScannr);
  void createScannrToMergedScannrMap(
    std::map<ScanId, ScanId>& scannrToMergedScannr);

  void writeClusterBins(unsigned int batchIdx,
    std::vector<pwiz::msdata::MSDataPtr>& msDataPtrs,
    std::vector<pwiz::msdata::SpectrumListSimplePtr>& spectrumListPtrMap);
};

} /* namespace maracluster */

#endif /* MARACLUSTER_MSFILEMERGER_H_ */
