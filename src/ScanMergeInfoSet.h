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
 
#ifndef MARACLUSTER_SCANMERGEINFOSET_H_
#define MARACLUSTER_SCANMERGEINFOSET_H_

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <climits>

#include <boost/foreach.hpp>

#include "ScanId.h"

namespace maracluster {

class ScanMergeInfoSet {
	public:
	  ScanMergeInfoSet() : minScannr_(UINT_MAX, UINT_MAX) {}
	  
	  inline size_t size() { return scans.size(); }
		void push_back(ScanId si) { 
		  scans.push_back(si);
		  if (si < minScannr_) minScannr_ = si;
		}
		std::vector<ScanId> scans;
		ScanId mergedScanId;
		ScanId minScannr_;
		
		inline static bool lowerScannr(const ScanMergeInfoSet& a, 
		                              const ScanMergeInfoSet& b) { 
      return (a.minScannr_ < b.minScannr_); 
    }
};

std::ostream& operator<<(std::ostream& os, const ScanMergeInfoSet& sms);
/* map: peptide -> ScanMergeInfoSet */
typedef std::map<std::string, ScanMergeInfoSet> SpectraMergeMap;

} /* namespace maracluster */

#endif /* MARACLUSTER_SCANMERGEINFOSET_H_ */
