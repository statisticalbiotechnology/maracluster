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
 
#ifndef BATCH_GLOBALS_H
#define BATCH_GLOBALS_H

#include <string>
#include <fstream>
#include <iostream>
#include <ctime>

#include "MyException.h"

class BatchGlobals {
 public:
  static unsigned int VERB; 
  
  static bool fileExists(const std::string& fileName);
  static void reportProgress(time_t& startTime, clock_t& startClock,
    size_t currentIt, size_t totalIt);
};

#endif // BATCH_GLOBALS_H
