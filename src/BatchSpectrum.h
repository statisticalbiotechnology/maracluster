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
 
#ifndef MARACLUSTER_BATCHSPECTRUM_H_
#define MARACLUSTER_BATCHSPECTRUM_H_

#ifdef FINGERPRINT_FILTER
  #define BATCH_SPECTRUM_NUM_STORED_PEAKS 80
#else
  #define BATCH_SPECTRUM_NUM_STORED_PEAKS 40
#endif

#include "ScanId.h"

namespace maracluster {

struct BatchSpectrum {
  ScanId scannr;
  unsigned int charge;
  float precMz, retentionTime;
  short fragBins[BATCH_SPECTRUM_NUM_STORED_PEAKS];
};

} /* namespace maracluster */

#endif /* MARACLUSTER_BATCHSPECTRUM_H_ */
