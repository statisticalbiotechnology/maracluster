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
 
#include <iostream>
#include <string>
#include <vector>

#include "pwiz/data/msdata/MSData.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/common/cv.hpp"
#include "pwiz/data/vendor_readers/ExtendedReaderList.hpp"
#include "Option.h"
#include "SpectrumHandler.h"
#include "MassChargeCandidate.h"

std::string spectrumOutFN_ = "";
std::string filePathOrig_ = "";
std::string filePathMerged_ = "";
std::string scannrString_ = "";
std::string clusterFile_ = "";

bool parseOptions(int argc, char **argv) {
  std::ostringstream intro;
  intro << "Usage:" << endl;
  intro << "  Option 1: extractspec -l spec_list -o spec_out" << std::endl;
  intro << "    spec_list is a tab delimited file of the form:" << std::endl;
  intro << "      filepath <tab> scannr <tab> ..." << std::endl;
  intro << "  Option 2: extractspec -i spec_in -m spec_merged_in \\";
  intro << std::endl;
  intro << "              -o spec_out -s scannr1,...,scannrN" << std::endl;
  intro << "    spec_in is an ms2 format file with original spectra.";
  intro << std::endl;
  intro << "    spec_merged_in is an ms2 format file with merged spectra."; 
  intro << std::endl;
  intro << "  spec_out is where the merged spectra will be written.";
  intro << std::endl;
  intro << "  (ensure to have read and write access on the files).";
  intro << std::endl;
  
  // init
  CommandLineParser cmd(intro.str());
  cmd.defineOption("l",
      "clusterFile",
      "File with filepaths and scannrs, separated by tabs",
      "filename");
  cmd.defineOption("i",
      "specIn",
      "File readable by ProteoWizard (e.g. ms2,mzML) with the original spectra",
      "filename");
  cmd.defineOption("o",
      "specOut",
      "File where you want the merged spectra to be written (ms2 format)",
      "filename");
  cmd.defineOption("m",
      "specMergeIn",
      "File readable by ProteoWizard (e.g. ms2, mzML) with the merged spectra",
      "filename");
  cmd.defineOption("s",
      "scannrs",
      "Scannrs as comma seperated string, e.g. 11,12,13",
      "int_list");
  
  // finally parse and handle return codes (display help etc...)
  cmd.parseArgs(argc, argv);
  
  // now query the parsing results
  if (cmd.optionSet("l")) clusterFile_ = cmd.options["l"];
  if (cmd.optionSet("o")) spectrumOutFN_ = cmd.options["o"];
  if (cmd.optionSet("i")) filePathOrig_ = cmd.options["i"];
  if (cmd.optionSet("m")) filePathMerged_ = cmd.options["m"];
  if (cmd.optionSet("s")) scannrString_ = cmd.options["s"];
  
  // if there are arguments left...
  if (cmd.arguments.size() > 0) {
    std::cerr << "Error: too many arguments." << std::endl;
    std::cerr << "Invoke with -h option for help." << std::endl;
    return 0; // ...error
  } else {
    if (!cmd.optionSet("i")) {
      std::cerr << "Error: one of the inputs is missing." << std::endl;
      std::cerr << "Invoke with -h option for help." << std::endl;
      return 0;
    }
  } 
  return true;
}

/* extractraw -i data/be/103111-Yeast-2hr-01.ms2 */
int main(int argc, char* argv[]) {
  try {
    if (parseOptions(argc, argv)) {
      pwiz::msdata::ExtendedReaderList readerList;
      
      pwiz::msdata::MSDataFile msd(filePathOrig_, &readerList);
      pwiz::msdata::SpectrumListPtr sl = msd.run.spectrumListPtr;

      for (unsigned int j = 0; j < sl->size(); ++j) {
        pwiz::msdata::SpectrumPtr s = sl->spectrum(j, true);
        if (s->cvParam(pwiz::cv::MS_ms_level).valueAs<int>() == 2) {
          std::vector<MassChargeCandidate> mccs;
          SpectrumHandler::getMassChargeCandidates(s, mccs);
          unsigned int scannr = SpectrumHandler::getScannr(s);
          std::cerr << "scannr: " << scannr << ", prec m/z: " << mccs.at(0).precMz << ", charge: " << mccs.at(0).charge << std::endl;
        }
      }
      return EXIT_SUCCESS;
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Caught unknown exception.\n";
  }

  return EXIT_FAILURE;
}

