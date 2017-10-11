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

#ifndef MARACLUSTER_PERCOLATORINTERFACE_H_
#define MARACLUSTER_PERCOLATORINTERFACE_H_

#include <cstdlib>
#include <stdexcept>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <vector>

#include <boost/foreach.hpp>

#include "ScanMergeInfo.h"
#include "SpectrumFileList.h"

namespace maracluster {

class PercolatorInterface {
 public:
  static unsigned int extractScannr(std::string id);
  static std::string extractCharge(std::string id);
  static int extractIntCharge(std::string id);
  
  static void parsePercOutfile(std::string percOutFN, SpectraMergeMap& peptideScanMap,
                      bool isDecoy, double qvalue_thresh);
  static void parsePercOutfile(std::string percOutFN, ScanPeptideMap& scanPeptideMap,
                      bool isDecoy);
  
  static void reducePercOutfile(SpectraMergeMap& peptideScanMap, 
    const std::string spectrumFileIn,
    SpectrumFileList& fileList, std::vector<ScanMergeInfoSet>& combineSets, 
    const std::string scanWeightsFN);
  static void reducePercInfile(std::string percInFN, std::string percInReducedFN, 
                      SpectraMergeMap& peptideScanMap);
};

} /* namespace maracluster */

#endif /* MARACLUSTER_PERCOLATORINTERFACE_H_ */
