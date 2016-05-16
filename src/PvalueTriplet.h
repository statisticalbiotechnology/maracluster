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
 
#ifndef PVALUE_TRIPLET_H
#define PVALUE_TRIPLET_H

#include "ScanId.h"

struct PvalueTriplet {
  ScanId scannr1, scannr2;
  float pval;
  
  PvalueTriplet() : scannr1(), scannr2(), pval(0.0) {}
  PvalueTriplet(ScanId _scannr1, ScanId _scannr2, float _pval) :
    scannr1(_scannr1), scannr2(_scannr2), pval(_pval) {}
  
  bool operator<(const PvalueTriplet& pt) const {
    return (pval < pt.pval) || (pval == pt.pval && scannr1 < pt.scannr1)
       || (pval == pt.pval && scannr1 == pt.scannr1 && scannr2 < pt.scannr2);
  }
  
  void readFromString(const char* f, char** next);
};

std::ostream& operator<<(std::ostream& stream, const PvalueTriplet& t);

#endif // PVALUE_TRIPLET_H
