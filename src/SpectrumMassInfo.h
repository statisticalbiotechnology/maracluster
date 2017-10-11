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
 
#ifndef MARACLUSTER_SPECTRUMMASSINFO_H_
#define MARACLUSTER_SPECTRUMMASSINFO_H_

#include <iostream>

namespace maracluster {

struct SpectrumMassInfo {
  SpectrumMassInfo(double _mass, unsigned int _spectrumIdx, unsigned int _charge, unsigned int _fileIdx = 0) : mass(_mass), spectrumIdx(_spectrumIdx), charge(_charge), fileIdx(_fileIdx) {}
  double mass;
  unsigned int spectrumIdx;
  unsigned int charge;
  unsigned int fileIdx;
};

// TODO: should we also sort by spectrum index to improve indexing?
inline static bool lessMass(const SpectrumMassInfo& a, const SpectrumMassInfo& b) { 
  return (a.mass < b.mass);
}

std::ostream& operator<<(std::ostream& os, const SpectrumMassInfo& mm) {
  os << std::setprecision(10) << "(" << mm.mass << "," << mm.spectrumIdx << "," << mm.fileIdx << "," << mm.charge << ")";
  return os;
}

} /* namespace maracluster */

#endif /* MARACLUSTER_SPECTRUMMASSINFO_H_ */
