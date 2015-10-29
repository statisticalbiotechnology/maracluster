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
 
#ifndef BATCH_PVALUES_H
#define BATCH_PVALUES_H

#include <iostream>
#include <vector>
#include <string>

#include "BatchGlobals.h"
#include "BinaryInterface.h"
#include "PvalueTriplet.h"

class BatchPvalues {
 public:
  BatchPvalues() : pvaluesFN_("") { }
  BatchPvalues(const std::string& pvaluesFN) : pvaluesFN_(pvaluesFN) { }
  
  void batchWrite(std::vector<PvalueTriplet>& pvalBuffer);
  
  inline std::string getPvaluesFN() const { return pvaluesFN_; }
 protected:
  std::string pvaluesFN_;
};

#endif // CASS_PVALUES_H
