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
 
#ifndef MARACLUSTER_MZINTENSITYPAIR_H_
#define MARACLUSTER_MZINTENSITYPAIR_H_

#include <vector>
#include <iostream>
#include <cstdint>

#include "pwiz/data/msdata/MSData.hpp"

namespace maracluster {

struct MZIntensityPair : public pwiz::msdata::MZIntensityPair {
  public:
    double multiplicity;
    
    MZIntensityPair(double mz_, double intensity_, double multiplicity_) {
      init(mz_, intensity_, multiplicity_);
    }
    
    MZIntensityPair(double mz_, double intensity_) {
      init(mz_, intensity_, 1.0);
    }
    
    MZIntensityPair(pwiz::msdata::MZIntensityPair mzi) {
      init(mzi.mz, mzi.intensity, 1.0);
    }
    
    MZIntensityPair() {
      init(0.0, 0.0, 1.0);
    }
  protected:
    void init(double mz_, double intensity_, double multiplicity_) {
      mz = mz_;
      intensity = intensity_;
      multiplicity = multiplicity_;
    }
};

std::ostream& operator<<(std::ostream& os, const MZIntensityPair& mzi);

} /* namespace maracluster */

#endif /* MARACLUSTER_MZINTENSITYPAIR_H_ */
