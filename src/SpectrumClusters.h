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
 
#ifndef MARACLUSTER_SPECTRUMCLUSTERS_H_
#define MARACLUSTER_SPECTRUMCLUSTERS_H_

#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include "Globals.h"
#include "SpectrumFiles.h"
#include "SpectrumFileList.h"
#include "PvalueTriplet.h"
#include "BinaryInterface.h"

namespace maracluster {

class SpectrumClusters {
 public:
  void printClusters(const std::vector<std::string>& pvalTreeFNs,
    const std::vector<double>& clusterThresholds, SpectrumFileList& fileList, 
    const std::string& scanInfoFN, const std::string& resultBaseFN,
    const std::string& scanTitleFN, bool addSpecIds);
  
  static std::string getClusterFN(const std::string resultBaseFN, double threshold);
  
 private:
  std::vector<ScanInfo> scanInfos_;
  std::map<ScanId, std::string> scanTitleMap_;
  
  void readPvalTree(const std::string& pvalTreeFN,
    std::vector<PvalueTriplet>& pvals);
  void readScanNrs(const std::string& scanInfoFN);
  
  void createClusterings(std::vector<PvalueTriplet>& pvals, 
    const std::vector<double>& clusterThresholds, SpectrumFileList& fileList,
    const std::string& resultBaseFN);
  void writeClusters(
    std::map<ScanId, std::vector<ScanId> >& clusters,
    SpectrumFileList& fileList, const std::string& resultFN);
  size_t writeSingletonClusters(std::set<ScanId>& seenScannrs,
    SpectrumFileList& fileList, std::ofstream& resultStream, size_t clusterIdx);
  void writeClusterSummary(
    std::vector<std::pair<size_t, size_t> >& clusterSizeCounts);
};

} /* namespace maracluster */

#endif /* MARACLUSTER_SPECTRUMCLUSTERS_H_ */
