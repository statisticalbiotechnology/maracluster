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
 
#ifndef MARACLUSTER_BATCHSPECTRUMFILES_H_
#define MARACLUSTER_BATCHSPECTRUMFILES_H_

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include "Globals.h"
#include "Pvalues.h"
#include "PvalueVectors.h"
#include "Spectrum.h"

#include "ScanId.h"
#include "PeakCounts.h"
#include "SpectrumFileList.h"
#include "SpectrumHandler.h"
#include "MSFileHandler.h"
#include "BinSpectra.h"
#include "BinaryInterface.h"
#include "MyException.h"

namespace maracluster {

struct ScanInfo {
  ScanInfo() : scanId(), minPrecMz(0.0), maxPrecMz(0.0) {}
  ScanInfo(const ScanId& si, const float precMz) :
    scanId(si), minPrecMz(precMz), maxPrecMz(precMz) {}
  
  ScanId scanId;
  float minPrecMz, maxPrecMz;
  
  bool operator<(const ScanInfo& si) const {
    return (scanId < si.scanId);
  }
  
  bool operator==(const ScanInfo& si) const {
    return (scanId == si.scanId);
  }
  
  bool operator!=(const ScanInfo& si) const {
    return (scanId != si.scanId);
  }
};

class SpectrumFiles {
 public:
  SpectrumFiles() : outputFolder_(""), chargeUncertainty_(0), datFolder_(""), addSpecIds_(false) {}
  SpectrumFiles(const std::string& precMzFileFolder,
                const std::string& datFolder) : 
      outputFolder_(precMzFileFolder), 
      datFolder_(datFolder), 
      chargeUncertainty_(0),
      addSpecIds_(false) {}
  SpectrumFiles(const std::string& precMzFileFolder, 
                const std::string& datFolder,
                const int chargeUncertainty,
                const bool addSpecIds) : 
      outputFolder_(precMzFileFolder),
      datFolder_(datFolder),
      chargeUncertainty_(chargeUncertainty),
      addSpecIds_(addSpecIds) {}
  
  void convertToDat(SpectrumFileList& fileList);
  void splitByPrecursorMz(SpectrumFileList& fileList,
      std::vector<std::string>& datFNs, const std::string& peakCountFN,
      const std::string& scanInfoFN, double precursorTolerance, 
      bool precursorToleranceDa);
  void splitByPrecursorMz(SpectrumFileList& fileList,
      const std::string& datFNFile, const std::string& peakCountFN,
      const std::string& scanInfoFN, double precursorTolerance, 
      bool precursorToleranceDa);
  
  void getPeakCountsAndPrecursorMzs(SpectrumFileList& fileList,
    std::vector<double>& precMzsAccumulated, const std::string& peakCountFN);
  
  void getBatchSpectra(const std::string& spectrumFN, 
    SpectrumFileList& fileList, std::vector<Spectrum>& localSpectra,
    std::vector<ScanInfo>& localScanInfos,
    std::vector<ScanIdExtended>& scanTitles);
  
  void getDatFNs(std::vector<double>& limits, 
    std::vector<std::string>& datFNs);
  static void writeDatFNsToFile(std::vector<std::string>& datFNs,
    const std::string& datFNFile);
  static void readDatFNsFromFile(const std::string& datFNFile,
    std::vector<std::string>& datFNs);
  
  static void writeScanTitlesToFile(std::vector<ScanIdExtended>& scanTitles,
    const std::string& scanTitlesFile);
    
  static void readPrecMzLimits(const std::string& scanInfoFN,
    std::map<ScanId, std::pair<float, float> >& precMzLimits);
  
  static bool limitsUnitTest();
  
 protected:
  std::string outputFolder_;
  std::string datFolder_;
  int chargeUncertainty_;
  bool addSpecIds_;
  
  virtual void getMassChargeCandidates(pwiz::msdata::SpectrumPtr s, 
    std::vector<MassChargeCandidate>& mccs, ScanId scanId);
  
  void writePrecMzs(const std::vector<double>& precMzs);
  void readPrecMzs(const std::string& precMzFN,
                             std::vector<double>& precMzs);
  
  void getPrecMzLimits(std::vector<double>& precMzs, 
    std::vector<double>& limits, double precursorTolerance, 
    bool precursorToleranceDa);
  int getPrecMzBin(double precMz, std::vector<double>& limits);
  
  void writeSplittedPrecursorMzFiles(SpectrumFileList& fileList, 
    std::vector<double>& limits, std::vector<std::string>& datFNs,
    const std::string& scanInfoFN);
  
  void appendBatchSpectra(
    std::vector< std::vector<Spectrum> >& batchSpectra,
    std::vector<std::string>& datFNs);

  void writePeakCounts(PeakCounts& peakCountsAccumulated, const std::string& peakCountFN);
  
  void convertAndWriteDatFiles(
    SpectrumFileList& fileList,
    const std::string& spectrumFN);
  void loadDatFiles(
    SpectrumFileList& fileList,
    const std::string& spectrumFN,
    std::vector<Spectrum>& localSpectra,
    std::vector<ScanInfo>& scanInfos);
  
  static std::string getFilename(const std::string& filepath);
  static std::string getDirectory(const std::string& filepath);
  static std::string getOutputFile(const std::string& filepath, 
                                   const std::string& outputFolder, 
                                   const std::string& newExtension);
  static void createDirectory(const boost::filesystem::path& dirPath);

};

} /* namespace maracluster */

#endif /* MARACLUSTER_BATCHSPECTRUMFILES_H_ */
