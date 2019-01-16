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

#ifndef MARACLUSTER_MARACLUSTER_H_
#define MARACLUSTER_MARACLUSTER_H_

#include <cstdio>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/filesystem.hpp>

#include "Option.h"
#include "PvalueCalculator.h"
#include "SpectrumFileList.h"

#include "Globals.h"
#include "Version.h"
#include "SpectrumFiles.h"
#include "SpectrumClusters.h"
#include "Spectra.h"
#include "PvalueVectors.h"
#include "Pvalues.h"

#include "MSFileExtractor.h"
#include "MSFileMerger.h"
#include "MSClusterMerge.h"

#include "PvalueFilterAndSort.h"
#include "SparseClustering.h"

namespace maracluster {

enum Mode { NONE, BATCH, PVALUE, UNIT_TEST, INDEX, CLUSTER, CONSENSUS, SEARCH, PROFILE_CONSENSUS, PROFILE_SEARCH };

class MaRaCluster {  
 public:
  MaRaCluster();
  ~MaRaCluster();
  
  bool parseOptions(int argc, char **argv);
  
  int run();
  
  virtual int mergeSpectra();
  
 protected:
  std::string greeter();
  std::string extendedGreeter(time_t& startTime);
  
  virtual int createIndex();
  int doClustering(const std::vector<std::string> pvalFNs, 
    std::vector<std::string> pvalTreeFNs, SpectrumFileList& fileList);
  
  Mode mode_;
  std::string call_;
  std::string percOutFN_;
  std::string fnPrefix_;
  std::string peakCountFN_;
  std::string datFNFile_;
  std::string scanInfoFN_;
  std::string pvaluesFN_;
  std::string clusterFileFN_;
  std::string pvalVecInFileFN_;
  std::string pvalueVectorsBaseFN_;
  std::string overlapBatchFileFN_;
  boost::filesystem::path outputPath;
  std::string outputFolder_;
  std::string spectrumBatchFileFN_;
  std::string spectrumInFN_;
  std::string spectrumOutFN_;
  std::string spectrumLibraryFN_;

  std::string matrixFN_;
  std::string resultTreeFN_;
  bool skipFilterAndSort_;
  bool writeAll_;
  std::vector<double> clusterThresholds_;
  double precursorTolerance_;
  bool precursorToleranceDa_;
  double dbPvalThreshold_; // logPval
  int chargeUncertainty_;
  size_t minConsensusClusterSize_;
};

} /* namespace maracluster */

#endif /* MARACLUSTER_MARACLUSTER_H_ */
