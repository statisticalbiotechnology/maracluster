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
 
#ifndef MARACLUSTER_BATCHPVALUES_H_
#define MARACLUSTER_BATCHPVALUES_H_

#include <iostream>
#include <vector>
#include <string>

#include "Globals.h"
#include "BinaryInterface.h"
#include "PvalueTriplet.h"

namespace maracluster {

class Pvalues {
 public:
  Pvalues() : pvaluesFN_("") { }
  Pvalues(const std::string& pvaluesFN) : pvaluesFN_(pvaluesFN) { }
  
  void batchWrite(std::vector<PvalueTriplet>& pvalBuffer,
                  const std::string& suffix = "");
  
  inline std::string getPvaluesFN() const { return pvaluesFN_; }
 protected:
  std::string pvaluesFN_;
};

} /* namespace maracluster */

#endif /* MARACLUSTER_BATCHPVALUES_H_ */
