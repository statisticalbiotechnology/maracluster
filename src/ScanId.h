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
 
#ifndef MARACLUSTER_SCANID_H_
#define MARACLUSTER_SCANID_H_

#include <boost/functional/hash.hpp>

namespace maracluster {

struct ScanId {
  unsigned int fileIdx, scannr;
  ScanId() : fileIdx(0u), scannr(0u) { }
  ScanId(unsigned int f, unsigned int s) : fileIdx(f), scannr(s) { }
  
  bool operator<(const ScanId& si) const {
    return (fileIdx < si.fileIdx) || 
           (fileIdx == si.fileIdx && scannr < si.scannr);
  }
  
  bool operator==(const ScanId& si) const {
    return (fileIdx == si.fileIdx && scannr == si.scannr);
  }
  
  bool operator!=(const ScanId& si) const {
    return (fileIdx != si.fileIdx || scannr != si.scannr);
  }
  
  void readFromString(const char* str, char** endptr);
};

std::ostream& operator<<(std::ostream& stream, const ScanId& si);
std::size_t hash_value(ScanId const& si);

} /* namespace maracluster */

#endif /* MARACLUSTER_SCANID_H_ */
