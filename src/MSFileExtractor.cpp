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
 
#include "MSFileExtractor.h"

namespace maracluster {

using pwiz::msdata::MSData;
using pwiz::msdata::MSDataFile;
using pwiz::msdata::SpectrumListSimplePtr;
using pwiz::msdata::SpectrumListSimple;
using pwiz::msdata::SpectrumListPtr;
using pwiz::msdata::SpectrumPtr;
using pwiz::msdata::Spectrum;
using pwiz::msdata::SelectedIon;

void MSFileExtractor::parseClusterFileForExtract(const std::string& clusterFile) {
  std::ifstream specListFN(clusterFile.c_str());
  std::string line;
  if (specListFN.is_open()) {
    combineSets_.push_back(ScanMergeInfoSet());
    while (getline(specListFN, line)) {
      std::string filepath;
      unsigned int scannr;
      
      std::istringstream iss(line);
      iss >> filepath >> scannr;
      if (iss.good()) {
        ScanId globalIdx = fileList_.getScanId(filepath, scannr);
        combineSets_.back().push_back(ScanMergeInfo(globalIdx));
      }
    }
  }
}

void MSFileExtractor::parseScannrStringForExtract(
    const std::string& scannrString, const std::string& filePathOrig, 
    const std::string& filePathMerged) {
  fileList_.addFile(filePathOrig);
  if (filePathMerged.size() > 0)
    fileList_.addFile(filePathMerged);
  
  combineSets_.push_back(ScanMergeInfoSet());
  
  std::stringstream ss(scannrString.c_str());
  std::string token;
  while (std::getline(ss, token, ',')) {
    unsigned int scannr = atoi(token.c_str());
    ScanId tmp = fileList_.getScanId(filePathOrig, scannr);
    combineSets_.back().push_back(ScanMergeInfo(tmp));
  }
}

void MSFileExtractor::getScanIdsByFile(
    std::vector< std::vector<ScanId> >& scanIdsByFile) {
  std::sort(combineSets_.front().scans.begin(), 
            combineSets_.front().scans.end(), ScanMergeInfo::lowerScannr);
  
  scanIdsByFile.resize(fileList_.getFilePaths().size());        
  BOOST_FOREACH (const ScanMergeInfo& scanMergeInfo, 
                 combineSets_.front().scans) {
    unsigned int fileIdx = fileList_.getFileIdx(scanMergeInfo.scannr);
    scanIdsByFile[fileIdx].push_back(scanMergeInfo.scannr);
  }
}

void MSFileExtractor::extractSpectra() {
  std::cerr << "Extracting spectra!\n";
  
  MSData msdExtracted;
  msdExtracted.id = msdExtracted.run.id = "extracted_spectra";
  SpectrumListSimplePtr extractedSpectra(new SpectrumListSimple);
  msdExtracted.run.spectrumListPtr = extractedSpectra;
  
  std::vector< std::vector<ScanId> > scanIdsByFile;
  getScanIdsByFile(scanIdsByFile);
  
  size_t idx = 0;
  for (size_t fileIdx = 0; fileIdx < scanIdsByFile.size(); ++fileIdx) {
    std::string filepath = fileList_.getFilePath(fileIdx);
    if ((fileIdx+1) % 20 == 0) {
      std::cerr << "Extracting spectra from file " << fileIdx+1 << "/" 
                << scanIdsByFile.size() << std::endl; 
    }
    MSReaderList readerList;
    MSDataFile msd(filepath, &readerList);
    SpectrumListPtr sl = msd.run.spectrumListPtr;
    
    for (size_t j = 0; j < sl->size(); ++j) {
      SpectrumPtr s = sl->spectrum(j, false);
      unsigned int scannr = SpectrumHandler::getScannr(s);
      ScanId scanId(fileIdx, scannr);
      if (std::find(scanIdsByFile[fileIdx].begin(), scanIdsByFile[fileIdx].end(), scanId) != scanIdsByFile[fileIdx].end()) {
        s = sl->spectrum(j, true);
        SpectrumHandler::setScannr(s, scanId);
        s->index = idx++;
        extractedSpectra->spectra.push_back(s);
      }
    }
  }
  
  writeMSData(msdExtracted, spectrumOutFN_); 
}

void MSFileExtractor::extractToBatchSpectrumList(
    std::vector<BatchSpectrum>& batchSpectra) {
  std::cerr << "Extracting spectra!\n";
  
  std::vector< std::vector<ScanId> > scanIdsByFile;
  getScanIdsByFile(scanIdsByFile);
  
  SpectrumListSimplePtr extractedSpectra(new SpectrumListSimple);
  for (size_t fileIdx = 0; fileIdx < scanIdsByFile.size(); ++fileIdx) {
    std::string filepath = fileList_.getFilePath(fileIdx);
    if ((fileIdx+1) % 20 == 0) {
      std::cerr << "Extracting spectra from file " << fileIdx+1 << "/" 
                << scanIdsByFile.size() << std::endl; 
    }
    MSReaderList readerList;
    MSDataFile msd(filepath, &readerList);
    SpectrumListPtr sl = msd.run.spectrumListPtr;
    
    BOOST_FOREACH (const ScanId scanId, scanIdsByFile[fileIdx]) {
      size_t result = getSpectrumIdxFromScannr(sl, scanId.scannr);
      SpectrumPtr s = sl->spectrum(result, true);
      
      std::vector<MZIntensityPair> mziPairs;
      SpectrumHandler::getMZIntensityPairs(s, mziPairs); 
      
      double retentionTime = SpectrumHandler::getRetentionTime(s);
      
      std::vector<MassChargeCandidate> mccs;
      int chargeUncertainty = 0;
      SpectrumHandler::getMassChargeCandidates(s, mccs, chargeUncertainty);
      
      BOOST_FOREACH (MassChargeCandidate& mcc, mccs) {
        int minCharge = (std::max)(static_cast<int>(mcc.charge) - static_cast<int>(chargeErrorTolerance_), 1);
        int maxCharge = mcc.charge + chargeErrorTolerance_;
        float precMz = SpectrumHandler::calcPrecMz(mcc.mass, mcc.charge);
        for (int charge = minCharge; charge <= maxCharge; ++charge) {
          BatchSpectrum bs;
          bs.precMz = precMz;
          bs.retentionTime = retentionTime;
          bs.charge = charge;
          bs.scannr = scanId;
          
          float precMass = SpectrumHandler::calcMass(precMz, charge);
          unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(precMass);
          std::vector<unsigned int> peakBins;
          BinSpectra::binBinaryTruncated(mziPairs, peakBins, numScoringPeaks, precMass);
          if (peakBins.size() >= PvalueCalculator::getMinScoringPeaks(precMass)) {
            std::copy(peakBins.begin(), peakBins.end(), bs.fragBins);
            batchSpectra.push_back(bs);
          }
        }
      }
    }
  }
}

} /* namespace maracluster */
