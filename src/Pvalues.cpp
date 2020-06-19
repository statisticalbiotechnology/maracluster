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
 
#include "Pvalues.h"

namespace maracluster {

void Pvalues::batchWrite(std::vector<PvalueTriplet>& pvalBuffer, 
    const std::string& suffix) {
  if (Globals::VERB > 4) {
    std::cerr << "Writing " << pvalBuffer.size() << " pvalues." << std::endl;
  }
  
  if (pvalBuffer.size() > 0)
#pragma omp critical (batch_write_pval)
  {  
    bool append = true;
    BinaryInterface::write<PvalueTriplet>(pvalBuffer, pvaluesFN_ + suffix, append);
  }
}

} /* namespace maracluster */
