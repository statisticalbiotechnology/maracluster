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
 
#ifndef MARACLUSTER_MSFILEEXTRACTOR_H_
#define MARACLUSTER_MSFILEEXTRACTOR_H_

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <cstdio>
#include <cstdint>

#include "pwiz/data/msdata/MSDataFile.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include "MSFileHandler.h"

namespace maracluster {

class MSFileExtractor : public MSFileHandler {
 public:  
  MSFileExtractor(std::string& spectrumOutFN) : 
    MSFileHandler(spectrumOutFN) {}
  
  void parseClusterFileForExtract(const std::string& clusterFile);
  
  void parseScannrStringForExtract(const std::string& scannrString,
    const std::string& filePathOrig, const std::string& filePathMerged);
  
  void extractSpectra();  
  void extractToBatchSpectrumList(std::vector<Spectrum>& batchSpectra);
  
 protected:  
  void getScanIdsByFile(std::vector< std::vector<ScanId> >& scanIdsByFile);
};

} /* namespace maracluster */

#endif /* MARACLUSTER_MSFILEEXTRACTOR_H_ */
