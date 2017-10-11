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

#include <cstdlib>

#include "MaRaCluster.h"

int main(int argc, char** argv) {
  maracluster::MaRaCluster maracluster;
  int retVal = EXIT_FAILURE; 
  
  try {
    if (maracluster.parseOptions(argc, argv)) {
    	retVal = maracluster.run();
	  }
  } catch (const std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << std::endl;
    retVal = EXIT_FAILURE;
  } catch(...) {
    std::cerr << "Unknown exception, contact the developer.." << std::endl;
    retVal = EXIT_FAILURE;
  }
  
  return retVal;
}
