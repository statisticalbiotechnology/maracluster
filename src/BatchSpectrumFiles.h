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

#include "Globals.h"
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

struct ScanInfo {
  ScanInfo() : scanId(), minPrecMz(0.0), maxPrecMz(0.0) {}
  ScanInfo(const ScanId& si, const float precMz) :
    scanId(si), minPrecMz(precMz), maxPrecMz(precMz) {}
  
  ScanId scanId;
  float minPrecMz, maxPrecMz;
};

class BatchSpectrumFiles {
 public:
  BatchSpectrumFiles() : precMzFileFolder_(""), chargeUncertainty_(0) {}
  BatchSpectrumFiles(const std::string& precMzFileFolder) : 
      precMzFileFolder_(precMzFileFolder), chargeUncertainty_(0) {}
  BatchSpectrumFiles(const std::string& precMzFileFolder, 
                     const int chargeUncertainty) : 
      precMzFileFolder_(precMzFileFolder), 
      chargeUncertainty_(chargeUncertainty) {}
  
  void splitByPrecursorMz(SpectrumFileList& fileList,
      std::vector<std::string>& datFNs, const std::string& peakCountFN,
      const std::string& scanInfoFN, double precursorTolerance, 
      bool precursorToleranceDa);
  void splitByPrecursorMz(SpectrumFileList& fileList,
      const std::string& datFNFile, const std::string& peakCountFN,
      const std::string& scanInfoFN, double precursorTolerance, 
      bool precursorToleranceDa);
  
  void writeDatFNsToFile(std::vector<std::string>& datFNs,
    const std::string& datFNFile);
  void getDatFNs(std::vector<double>& limits, std::vector<std::string>& datFNs);
  void readDatFNsFromFile(const std::string& datFNFile,
    std::vector<std::string>& datFNs);
  void readPrecMzLimits(const std::string& scanInfoFN,
    std::map<ScanId, std::pair<float, float> >& precMzLimits);
  
  void getBatchSpectra(const std::string& spectrumFN, 
    SpectrumFileList& fileList, std::vector<BatchSpectrum>& localSpectra,
    std::vector<ScanInfo>& localScanInfos);
  
  static bool limitsUnitTest();
  
 protected:
  std::string precMzFileFolder_;
  int chargeUncertainty_;
  
  void writePrecMzs(const std::vector<double>& PrecMzs);
  void readPrecMzs(const std::string& precMzFN,
                             std::vector<double>& PrecMzs);
  
  void getPrecMzLimits(std::vector<double>& precMzs, 
    std::vector<double>& limits, double precursorTolerance, 
    bool precursorToleranceDa);
  int getPrecMzBin(double precMz, std::vector<double>& limits);
  
  void getPeakCountsAndPrecursorMzs(SpectrumFileList& fileList,
    std::vector<double>& precMzsAccumulated, const std::string& peakCountFN);
  void writeSplittedPrecursorMzFiles(SpectrumFileList& fileList, 
    std::vector<double>& limits, std::vector<std::string>& datFNs,
    const std::string& scanInfoFN);
  
  void appendBatchSpectra(
    std::vector< std::vector<BatchSpectrum> >& batchSpectra,
    std::vector<std::string>& datFNs);
  void writePeakCounts(PeakCounts& peakCountsAccumulated, const std::string& peakCountFN);
  
  static std::string getFilename(const std::string& filepath);
  static std::string getDirectory(const std::string& filepath);
  
};

#endif // CASS_SPECTRUM_FILES_H
