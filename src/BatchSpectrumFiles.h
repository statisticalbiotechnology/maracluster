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
 
#ifndef CASS_SPECTRUM_FILES_H
#define CASS_SPECTRUM_FILES_H

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include <boost/foreach.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/data/msdata/MSDataMerger.hpp"

#include "BatchGlobals.h"
#include "BatchStatement.h"
#include "BatchPvalues.h"
#include "BatchPvalueVectors.h"
#include "BatchSpectrum.h"

#include "PeakCounts.h"
#include "SpectrumFileList.h"
#include "SpectrumHandler.h"
#include "MSFileHandler.h"
#include "BinSpectra.h"
#include "BinaryInterface.h"

class BatchSpectrumFiles {
 public:
  BatchSpectrumFiles() : precMassFileFolder_(""), chargeUncertainty_(0) {}
  BatchSpectrumFiles(const std::string& precMassFileFolder) : 
      precMassFileFolder_(precMassFileFolder), chargeUncertainty_(0) {}
  BatchSpectrumFiles(const std::string& precMassFileFolder, 
                     const int chargeUncertainty) : 
      precMassFileFolder_(precMassFileFolder), 
      chargeUncertainty_(chargeUncertainty) {}
  
  void splitByPrecursorMass(SpectrumFileList& fileList,
      std::vector<std::string>& datFNs, const std::string& peakCountFN,
      const std::string& scanNrsFN, double precursorTolerance, 
      bool precursorToleranceDa);
  void splitByPrecursorMass(SpectrumFileList& fileList,
      const std::string& datFNFile, const std::string& peakCountFN,
      const std::string& scanNrsFN, double precursorTolerance, 
      bool precursorToleranceDa);
  
  void writeDatFNsToFile(std::vector<std::string>& datFNs,
    const std::string& datFNFile);
  void getDatFNs(std::vector<double>& limits, std::vector<std::string>& datFNs);
  void readDatFNsFromFile(const std::string& datFNFile,
    std::vector<std::string>& datFNs);
  
  void writeScannrs(SpectrumFileList& fileList, 
                    const std::string& scanNrsFN);
  
  static bool limitsUnitTest();
  
 protected:
  std::string precMassFileFolder_;
  int chargeUncertainty_;
  
  void writePrecMasses(const std::vector<double>& precMasses);
  void readPrecMasses(const std::string& precMassFN,
                             std::vector<double>& precMasses);
  
  void getPrecMassLimits(std::vector<double>& precMasses, 
    std::vector<double>& limits, double precursorTolerance, 
    bool precursorToleranceDa);
  int getPrecMassBin(double precMass, std::vector<double>& limits);
  
  void getPeakCountsAndPrecursorMasses(SpectrumFileList& fileList,
    std::vector<double>& precMassesAccumulated, const std::string& peakCountFN);
  void writeSplittedPrecursorMassFiles(SpectrumFileList& fileList, 
    std::vector<double>& limits, std::vector<std::string>& datFNs,
    const std::string& scanNrsFN);
  
  void appendBatchSpectra(
    std::vector< std::vector<BatchSpectrum> >& batchSpectra,
    std::vector<std::string>& datFNs);
  void writePeakCounts(PeakCounts& peakCountsAccumulated, const std::string& peakCountFN);
  
  static std::string getFilename(const std::string& filepath);
  static std::string getDirectory(const std::string& filepath);
  
};

#endif // CASS_SPECTRUM_FILES_H
