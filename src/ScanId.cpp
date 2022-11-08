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

#include "ScanId.h"

namespace maracluster {

void ScanId::readFromString(const char* f, char** next) {
  fileIdx = strtoul(f, next, 0); f = *next;
  scannr = strtoul(f, next, 0);
}

std::ostream& operator<<(std::ostream& stream, const ScanId& si) {
  stream << si.fileIdx << " " << si.scannr;
  return stream;
}

std::size_t hash_value(ScanId const& si) {
  boost::hash<int> hasher;
  return hasher((si.fileIdx+1) * 1000000 + si.scannr);
}

std::ostream& operator<<(std::ostream& stream, const ScanIdExtended& si) {
  stream << si.scanId.fileIdx << " " << si.scanId.scannr << " " << si.title;
  return stream;
}

} /* namespace maracluster */
