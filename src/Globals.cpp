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
 
#include "Globals.h"
#include "MyException.h"

namespace maracluster {

unsigned int Globals::VERB = 3;

bool Globals::fileExists(const std::string& fileName) {
  std::ifstream infile(fileName.c_str());
  return infile.good();
}

bool Globals::fileIsEmpty(const std::string& fileName) {
  std::ifstream in(fileName.c_str(), std::ios::ate | std::ios::binary);
  if (in.is_open()) {
    return static_cast<long long>(in.tellg()) == 0ll;
  } else {
    return true;
  }
}

std::string Globals::getDirectory(const std::string& filepath) {
  unsigned found = filepath.find_last_of("/\\");
  return filepath.substr(0,found);
}

std::string Globals::getFilename(const std::string& filepath) {
  unsigned found = filepath.find_last_of("/\\");
  return filepath.substr(found+1);
}

std::string Globals::getOutputFile(const std::string& filepath, 
    const std::string& outputFolder, const std::string& newExtension) {
  std::string filename = getFilename(filepath);
  unsigned found = filename.find_last_of(".");
  return outputFolder + "/" + filename.substr(0,found) + newExtension;
}

void Globals::createDirectory(const boost::filesystem::path& dirPath) {
  if (boost::filesystem::exists(dirPath)) return;
  
  boost::system::error_code returnedError;
  boost::filesystem::create_directories(dirPath, returnedError );
  if (!boost::filesystem::exists(dirPath)) {
    std::stringstream ss;
    ss << "(SpectrumFiles.cpp) error creating folder " << dirPath 
       << " (" << returnedError.message() << ")" << std::endl;
    throw MyException(ss);
  }
}

void Globals::reportProgress(time_t& startTime, clock_t& startClock,
    size_t currentIt, size_t totalIt) {
  time_t elapsedTime;
  time(&elapsedTime);
  clock_t elapsedClock = clock();
  
  double diff = difftime(elapsedTime, startTime);
  
  unsigned int timeElapsedMin = static_cast<unsigned int>(diff/60);
  unsigned int timeElapsedSecMod = 
      static_cast<unsigned int>(diff - timeElapsedMin * 60);
  
  double elapsedCpuTime = (elapsedClock - startClock) / (double)CLOCKS_PER_SEC;
  std::cerr << "  Elapsed time: " << elapsedCpuTime << " cpu seconds " <<
               "or " << timeElapsedMin << " min " << timeElapsedSecMod << 
               " sec wall time." << std::endl;
  
  double timeLeftSec = (diff / (currentIt+1)) * (totalIt - (currentIt+1));
  unsigned int timeLeftMin = static_cast<unsigned int>(timeLeftSec/60);
  unsigned int timeLeftSecMod = 
      static_cast<unsigned int>(timeLeftSec - timeLeftMin * 60);
  std::cerr << "  Estimated time remaining: " << timeLeftMin << " min " <<
               timeLeftSecMod << " sec wall time." << std::endl;
}


} /* namespace maracluster */
