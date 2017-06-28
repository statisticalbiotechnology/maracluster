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
 
#include "SpectrumFileList.h"
#include <climits>

ScanId SpectrumFileList::getScanId(const std::string& filePath, 
    const unsigned int scannr) {
  addFile(filePath);
  return ScanId(fileIndexMap_[filePath], scannr);
}

void SpectrumFileList::initFromFile(const std::string& fileListFN) {
  std::ifstream specListFN(fileListFN.c_str());
  std::string line;
  if (specListFN.is_open()) {
    while (getline(specListFN, line)) {
      std::string filePath;
      
      std::istringstream iss(line);
      iss >> filePath;
      addFile(filePath);
    }
    
    if (fileIndexVector_.size() == 0) {
      std::stringstream ss;
      ss << "(SpectrumFileList.cpp) no input files detected in " << fileListFN << std::endl;
      throw MyException(ss);
    }
  } else {
    std::stringstream ss;
    ss << "(SpectrumFileList.cpp) could not open " << fileListFN << std::endl;
    throw MyException(ss);
  }
}

void SpectrumFileList::addFile(const std::string& filePath) {
  if (fileIndexMap_.find(filePath) == fileIndexMap_.end()) {
    fileIndexMap_[filePath] = fileIndexVector_.size();
    fileIndexVector_.push_back(filePath);
  }
}
