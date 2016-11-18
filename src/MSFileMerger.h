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
 
#ifndef MS_FILE_MERGER_H
#define MS_FILE_MERGER_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <cstdio>

#include "pwiz/data/msdata/MSDataFile.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include "Globals.h"
#include "MSFileHandler.h"
#include "MSClusterMerge.h"
#include "ClusterMerge.h"
#include "InterpolationMerge.h"
#include "RankMerge.h"
#include "ScanId.h"

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
  
  void parseClusterFileForSingleFileMerge(
      const std::string& clusterFN, const std::string& spectrumInFN,
      const std::string& scanWeightsFN);
  
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
  std::vector<pwiz::msdata::SourceFilePtr> sourceFilePtrs_;
  
  static std::string getPartFN(const std::string& outputFN, 
                               const std::string& partString);
  
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
    pwiz::msdata::SpectrumListSimplePtr mergedSpectra);
    
  void mergeSplitSpecFiles();
  
  void splitSpecFilesByConsensusSpec(
    std::map<ScanId, ScanId>& scannrToMergedScannr);
  void createScannrToMergedScannrMap(
    std::map<ScanId, ScanId>& scannrToMergedScannr);
  
  void writeClusterBins(unsigned int batchIdx,
    std::vector<pwiz::msdata::SpectrumListSimplePtr>& spectrumListPtrMap);
};

#endif // MS_FILE_MERGER_H
