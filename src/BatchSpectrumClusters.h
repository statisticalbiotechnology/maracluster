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
 
#ifndef MARACLUSTER_BATCHSPECTRUMCLUSTERS_H_
#define MARACLUSTER_BATCHSPECTRUMCLUSTERS_H_

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
#include "BatchSpectrumFiles.h"
#include "SpectrumFileList.h"
#include "PvalueTriplet.h"
#include "ScanMergeInfo.h"
#include "BinaryInterface.h"

namespace maracluster {

class BatchSpectrumClusters {
 public:
  void printClusters(const std::vector<std::string>& pvalTreeFNs,
    const std::vector<double>& clusterThresholds, SpectrumFileList& fileList, 
    const std::string& scanInfoFN, const std::string& scanDescFN,
    const std::string& resultBaseFN);
  
  static std::string getClusterFN(const std::string resultBaseFN, double threshold);
  
  static bool scanDescReadUnitTest();
 private:
  ScanPeptideMap scanPeptideMap_;
  
  void readPvalTree(const std::string& pvalTreeFN,
    std::vector<PvalueTriplet>& pvals);
  void createScanDescriptionMap(
    const std::string& scanInfoFN, const std::string& scanDescFN,
    SpectrumFileList& fileList);
  void readScanNrs(const std::string& scanInfoFN);
  void readScanDescs(const std::string& scanDescFN, SpectrumFileList& fileList);
  
  void createClusterings(std::vector<PvalueTriplet>& pvals, 
    const std::vector<double>& clusterThresholds, SpectrumFileList& fileList,
    const std::string& resultBaseFN);
  void writeClusters(
    std::map<ScanId, std::vector<ScanId> >& clusters,
    SpectrumFileList& fileList, const std::string& resultFN);
  size_t writeSingletonClusters(std::set<ScanId>& seenScannrs,
    SpectrumFileList& fileList, std::ofstream& resultStream);
  void writeClusterSummary(
    std::vector<std::pair<size_t, size_t> >& clusterSizeCounts);
};

} /* namespace maracluster */

#endif /* MARACLUSTER_BATCHSPECTRUMCLUSTERS_H_ */
