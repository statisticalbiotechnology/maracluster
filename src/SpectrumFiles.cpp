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
 
#include "SpectrumFiles.h"

namespace maracluster {

using pwiz::msdata::SpectrumListPtr;
using pwiz::msdata::MSDataFile;
using pwiz::msdata::SpectrumPtr;

void SpectrumFiles::convertToDat(
    SpectrumFileList& fileList) {
  if (Globals::VERB > 1) {
    std::cerr << "Converting input files to .dat binary format" << std::endl;
  }
  
  SpectrumFiles::createDirectory(datFolder_);
  
  std::vector<std::string> spectrumFNs = fileList.getFilePaths();
#pragma omp parallel for schedule(dynamic, 1)                
  for (int fileIdx = 0; fileIdx < static_cast<int>(spectrumFNs.size()); ++fileIdx) {
    std::string spectrumFN = spectrumFNs[fileIdx];
    if (Globals::VERB > 1) {
      std::cerr << "  Processing " << spectrumFN << 
          " (" << (fileIdx+1)*100/spectrumFNs.size() << "%)." << std::endl;
    }
    convertAndWriteDatFiles(fileList, spectrumFN);
  }
}

void SpectrumFiles::splitByPrecursorMz(
    SpectrumFileList& fileList, std::vector<std::string>& datFNs,
    const std::string& peakCountFN, const std::string& scanInfoFN,
    const std::string& scanTitleFN, double precursorTolerance, 
    bool precursorToleranceDa) {
  if (Globals::VERB > 1) {
    std::cerr << "Splitting spectra by precursor Mz" << std::endl;
  }
  
  std::vector<double> precMzsAccumulated;
  getPeakCountsAndPrecursorMzs(fileList, precMzsAccumulated, peakCountFN);
  
  //TsvInterface::write<double>(precMzsAccumulated, peakCountFN_, false);
  //TsvInterface::read<double>(peakCountFN_, precMzsAccumulated);
  
  std::vector<double> limits;
  getPrecMzLimits(precMzsAccumulated, limits, precursorTolerance, 
                    precursorToleranceDa);
  getDatFNs(limits, datFNs);
  writeSplittedPrecursorMzFiles(fileList, limits, datFNs, scanInfoFN, scanTitleFN);
}

void SpectrumFiles::splitByPrecursorMz(SpectrumFileList& fileList,
    const std::string& datFNFile, const std::string& peakCountFN,
    const std::string& scanInfoFN, const std::string& scanTitleFN, 
    double precursorTolerance, bool precursorToleranceDa) {  
  std::vector<std::string> datFNs;
  splitByPrecursorMz(fileList, datFNs, peakCountFN, scanInfoFN, scanTitleFN, 
                       precursorTolerance, precursorToleranceDa);
  
  bool append = false;
  TsvInterface::write<std::string>(datFNs, datFNFile, append);
}

void SpectrumFiles::getPeakCountsAndPrecursorMzs(
    SpectrumFileList& fileList, 
    std::vector<double>& precMzsAccumulated,
    const std::string& peakCountFN) {
  if (Globals::VERB > 1) {
    std::cerr << "Accumulating peak counts and precursor Mzs" << std::endl;
  }
  
  PeakCounts peakCountsAccumulated;
  
  std::vector<std::string> spectrumFNs = fileList.getFilePaths();
#pragma omp parallel for schedule(dynamic, 1)  
  for (int fileIdx = 0; fileIdx < static_cast<int>(spectrumFNs.size()); ++fileIdx) {
    std::string spectrumFN = spectrumFNs[fileIdx];
    if (Globals::VERB > 1 && fileIdx % 100 == 0) {
      std::cerr << "  Processing file " << fileIdx+1 << "/" << spectrumFNs.size() << 
          " (" << (fileIdx+1)*100/spectrumFNs.size() << "%)." << std::endl;
    }

    PeakCounts peakCounts;
    std::vector<double> precMzs;
    
    std::vector<Spectrum> localSpectra;
    std::vector<ScanInfo> scanInfos;
    std::vector<ScanIdExtended> scanTitles;
    loadDatFiles(fileList, spectrumFN, localSpectra, scanInfos, scanTitles);
    
    unsigned int lastCharge = 0;
    ScanId lastScannr;
    BOOST_FOREACH (Spectrum& spectrum, localSpectra) {
      precMzs.push_back(spectrum.precMz);
      unsigned int charge = (std::min)(spectrum.charge, peakCounts.getMaxCharge());
      double mass = SpectrumHandler::calcMass(spectrum.precMz, spectrum.charge);
      
      // only count each spectrum once per charge state
      if (charge == lastCharge && spectrum.scannr == lastScannr) {
        continue;
      }
      
      unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(mass);
      std::vector<unsigned int> peakBins(spectrum.fragBins, spectrum.fragBins + SPECTRUM_NUM_STORED_PEAKS);
      peakCounts.addSpectrum(peakBins, spectrum.precMz, charge, mass, numScoringPeaks);
      
      lastCharge = charge;
      lastScannr = spectrum.scannr;
    }
  #pragma omp critical (add_to_peakcount)  
    {
      peakCountsAccumulated.add(peakCounts);
      precMzsAccumulated.insert( precMzsAccumulated.end(), precMzs.begin(), precMzs.end() );
    }
  }
  writePeakCounts(peakCountsAccumulated, peakCountFN);
  
  std::sort(precMzsAccumulated.begin(), precMzsAccumulated.end());
}

void SpectrumFiles::writeSplittedPrecursorMzFiles(
    SpectrumFileList& fileList, 
    std::vector<double>& limits,
    std::vector<std::string>& datFNs,
    const std::string& scanInfoFN,
    const std::string& scanTitleFN) {
  if (Globals::VERB > 1) {
    std::cerr << "Dividing spectra in " << limits.size() << 
                 " bins of ~2 CPU hours each." << std::endl;
  }
  
  std::vector<std::string> spectrumFNs = fileList.getFilePaths();
#pragma omp parallel for schedule(dynamic, 1)                
  for (int fileIdx = 0; fileIdx < static_cast<int>(spectrumFNs.size()); ++fileIdx) {
    std::string spectrumFN = spectrumFNs[fileIdx];
    if (Globals::VERB > 1 && fileIdx % 100 == 0) {
      std::cerr << "  Processing file " << fileIdx+1 << "/" << spectrumFNs.size() << 
          " (" << (fileIdx+1)*100/spectrumFNs.size() << "%)." << std::endl;
    }
    
    std::vector<Spectrum> localSpectra;
    std::vector<ScanInfo> scanInfos;
    std::vector<ScanIdExtended> scanTitles;
    loadDatFiles(fileList, spectrumFN, localSpectra, scanInfos, scanTitles);
    
    std::vector< std::vector<Spectrum> > batchSpectra(limits.size());
    BOOST_FOREACH (Spectrum& bs, localSpectra) {
      int precBin = getPrecMzBin(bs.precMz, limits);
      batchSpectra[precBin].push_back(bs);
    }
    
  #pragma omp critical (add_to_datfiles)
    {
      appendBatchSpectra(batchSpectra, datFNs);
    }
  #pragma omp critical (write_scannrs)
    {
      bool append = true;
      BinaryInterface::write<ScanInfo>(scanInfos, scanInfoFN, append);
      if (addSpecIds_) {
        TsvInterface::write<ScanIdExtended>(scanTitles, scanTitleFN, append);
      }
    }
  }
}


void SpectrumFiles::getBatchSpectra(
    const std::string& spectrumFN, SpectrumFileList& fileList,
    std::vector<Spectrum>& localSpectra, 
    std::vector<ScanInfo>& scanInfos,
    std::vector<ScanIdExtended>& scanTitles) {
  if ( !boost::filesystem::exists( spectrumFN ) ) {
    std::cerr << "Ignoring missing file " << spectrumFN << std::endl;
    return;
  }
   
  SpectrumListPtr specList;
  MSReaderList readerList;
  MSDataFile msd(spectrumFN, &readerList);
  specList = msd.run.spectrumListPtr;
  
  size_t numSpectra = specList->size();
  for (size_t i = 0; i < numSpectra; ++i) {
    SpectrumPtr s = specList->spectrum(i, true);
    if (!SpectrumHandler::isMs2Scan(s)) continue;
    
    std::vector<MZIntensityPair> mziPairs;
    SpectrumHandler::getMZIntensityPairs(s, mziPairs); 
    
    double retentionTime = SpectrumHandler::getRetentionTime(s);
    unsigned int scannr = SpectrumHandler::getScannr(s);
    ScanInfo scanInfo;
    scanInfo.scanId = fileList.getScanId(spectrumFN, scannr);
    
    if (addSpecIds_) {
      ScanIdExtended scanTitle;
      scanTitle.scanId = scanInfo.scanId;
      scanTitle.title = SpectrumHandler::getScanTitle(s);
      scanTitles.push_back(scanTitle);
    }
    
    std::vector<MassChargeCandidate> mccs;
    getMassChargeCandidates(s, mccs, scanInfo.scanId);
    
    BOOST_FOREACH (MassChargeCandidate& mcc, mccs) {
      for (int isotopeTolerance = 0; isotopeTolerance <= 0; ++isotopeTolerance) {
        double mass = mcc.mass + isotopeTolerance;
        std::vector<unsigned int> peakBins;
        unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(mass);
        BinSpectra::binBinaryTruncated(mziPairs, peakBins, 
          numScoringPeaks, mcc.mass);
        
        if (peakBins.size() >= PvalueCalculator::getMinScoringPeaks(mass)) {
          Spectrum bs;
          peakBins.resize(SPECTRUM_NUM_STORED_PEAKS, 0u);
          std::copy(peakBins.begin(), peakBins.end(), bs.fragBins);
          bs.precMz = static_cast<float>(SpectrumHandler::calcPrecMz(mass, mcc.charge));
          bs.retentionTime = static_cast<float>(retentionTime);
          bs.charge = mcc.charge;
          bs.scannr = scanInfo.scanId;
          
          if (scanInfo.minPrecMz == 0.0 || bs.precMz < scanInfo.minPrecMz)
            scanInfo.minPrecMz = bs.precMz;
          if (scanInfo.maxPrecMz == 0.0 || bs.precMz > scanInfo.maxPrecMz)
            scanInfo.maxPrecMz = bs.precMz;
          
          localSpectra.push_back(bs);
        }
      }
    }
    scanInfos.push_back(scanInfo);
  }
}

void SpectrumFiles::appendBatchSpectra(
    std::vector< std::vector<Spectrum> >& batchSpectra,
    std::vector<std::string>& datFNs) {
  bool append = true;
  for (size_t i = 0; i < datFNs.size(); ++i) {
    BinaryInterface::write<Spectrum>(batchSpectra[i], datFNs[i], append);
  }
}

void SpectrumFiles::getDatFNs(std::vector<double>& limits, 
    std::vector<std::string>& datFNs) {
  int lastLimit = -1;
  int counter = 0;
  BOOST_FOREACH (double limit, limits) {
    int intLimit = static_cast<int>(limit);
    std::string filepath = outputFolder_ + "/" + 
        boost::lexical_cast<std::string>(intLimit);
    if (intLimit == lastLimit) {
      filepath += "_" + boost::lexical_cast<std::string>(++counter);
    } else {
      lastLimit = intLimit;
      counter = 0;
    }
    filepath += ".dat";
    datFNs.push_back(filepath);
  }
}

void SpectrumFiles::loadDatFiles(
    SpectrumFileList& fileList,
    const std::string& spectrumFN, 
    std::vector<Spectrum>& localSpectra,
    std::vector<ScanInfo>& scanInfos,
    std::vector<ScanIdExtended>& scanTitles) {
  std::string datFile = SpectrumFiles::getOutputFile(spectrumFN, datFolder_, ".dat");
  std::string scanInfoFile = SpectrumFiles::getOutputFile(spectrumFN, datFolder_, ".scan_info.dat");
  std::string scanTitleFN = SpectrumFiles::getOutputFile(spectrumFN, datFolder_, ".scan_titles.txt");
  if (!boost::filesystem::exists(datFile) || !boost::filesystem::exists(scanInfoFile)) {
    std::stringstream ss;
    ss << "(SpectrumFiles.cpp) missing dat file " << datFile
       << " and/or " <<  scanInfoFile << std::endl;
    throw MyException(ss);
  }
  
  BinaryInterface::read<Spectrum>(datFile, localSpectra);
  BinaryInterface::read<ScanInfo>(scanInfoFile, scanInfos);
  if (addSpecIds_) {
    TsvInterface::read<ScanIdExtended>(scanTitleFN, scanTitles);
  }
  
  // we need to update the file index, because we cannot guarantee that the 
  // input files are in the same order or if new files were added in between!
  unsigned int fileIdx = fileList.getFileIdx(spectrumFN);
  BOOST_FOREACH (ScanInfo& si, scanInfos) {
    si.scanId.fileIdx = fileIdx;
  }
  BOOST_FOREACH (Spectrum& s, localSpectra) {
    s.scannr.fileIdx = fileIdx;
  }
  if (addSpecIds_) {
    BOOST_FOREACH (ScanIdExtended& s, scanTitles) {
      s.scanId.fileIdx = fileIdx;
    }
  }
}

void SpectrumFiles::convertAndWriteDatFiles(
    SpectrumFileList& fileList,
    const std::string& spectrumFN) {
  std::string datFile = SpectrumFiles::getOutputFile(spectrumFN, datFolder_, ".dat");
  std::string scanInfoFile = SpectrumFiles::getOutputFile(spectrumFN, datFolder_, ".scan_info.dat");
  std::string scanTitleFN = SpectrumFiles::getOutputFile(spectrumFN, datFolder_, ".scan_titles.txt");
  if (boost::filesystem::exists(datFile) && boost::filesystem::exists(scanInfoFile)) {
    return;
  }
  
  std::vector<Spectrum> localSpectra;
  std::vector<ScanInfo> scanInfos;
  std::vector<ScanIdExtended> scanTitles;
  getBatchSpectra(spectrumFN, fileList, localSpectra, scanInfos, scanTitles);
  
  bool append = false;
  BinaryInterface::write<Spectrum>(localSpectra, datFile, append);
  BinaryInterface::write<ScanInfo>(scanInfos, scanInfoFile, append);
  if (addSpecIds_) {
    TsvInterface::write<ScanIdExtended>(scanTitles, scanTitleFN, append);
  }
}

void SpectrumFiles::readPrecMzLimits(const std::string& scanInfoFN,
    std::map<ScanId, std::pair<float, float> >& precMzLimits) {
  if (Globals::VERB > 2) {
    std::cerr << "Reading precursor m/z limits." << std::endl;
  }
  std::vector<ScanInfo> scanInfos;
  BinaryInterface::read<ScanInfo>(scanInfoFN, scanInfos);
  
  BOOST_FOREACH (const ScanInfo& si, scanInfos) {
    precMzLimits[si.scanId] = std::make_pair(si.minPrecMz, si.maxPrecMz);
  }
}

void SpectrumFiles::getPrecMzLimits(std::vector<double>& precMzs, 
    std::vector<double>& limits, double precursorTolerance, 
    bool precursorToleranceDa) {
  float pvecCost = 0.7f; // computation time for 1 p-value vector in ms on 1 core
  float pvalCost = 0.001f; // computation time for 1 p-value pair in ms on 1 core
  //float maxCost = 0.1f*60.0f*1000.0f; // test for splitted precursor bins
  float maxCost = 120.0f*60.0f*1000.0f; // max computation time = 120 CPU minutes
  
  unsigned long long numComparisons = 0uL;
  if (precMzs.size() > 0) {
    size_t lowerBoundIdx = 0;
    limits.push_back(precMzs[0]);
    float curCost = 0.0f;
    for (size_t idx = 0; idx < precMzs.size(); ++idx) {
      curCost += pvecCost;
      double lowerBound = PvalueVectors::getLowerBound(precMzs[idx], 
          precursorTolerance, precursorToleranceDa);
      while (precMzs[lowerBoundIdx] < lowerBound) {
        ++lowerBoundIdx;
      }
      numComparisons += (idx - lowerBoundIdx);
      curCost += (idx - lowerBoundIdx)*pvalCost;
      // ensure that each bin contains at least two precursor tolerance widths 
      // such that the overlap regions do not overlap
      double lowerLowerBound = PvalueVectors::getLowerBound(lowerBound, 
          precursorTolerance, precursorToleranceDa);
      if (curCost > maxCost && lowerLowerBound > limits.back()) {
        curCost = 0.0f;
        limits.push_back(precMzs[idx]);
      }
    } 
  }
  
  if (Globals::VERB > 1) {
    std::cerr << "Estimated number of pair comparisons: " 
              << numComparisons << std::endl;
  }
}

int SpectrumFiles::getPrecMzBin(double precMz, std::vector<double>& limits) {
  int bin = std::upper_bound(limits.begin(), limits.end(), precMz) - limits.begin() - 1;
  return std::min(std::max(bin, 0), static_cast<int>(limits.size()) - 1);
}

void SpectrumFiles::writePeakCounts(PeakCounts& peakCountsAccumulated, 
                                         const std::string& peakCountFN) {
  if (Globals::VERB > 2) {
    std::cerr << "Writing peak counts to file" << std::endl;
  }
  std::string serializedPeakCounts;
  PeakCounts::serializePeakCounts(peakCountsAccumulated, serializedPeakCounts);
  
  std::ofstream peakCountStream(peakCountFN.c_str(), std::ios_base::binary | std::ios_base::out);
  if (peakCountStream.is_open()) {
    peakCountStream << serializedPeakCounts;
  }
  
  if (Globals::VERB > 2) {
    std::cerr << "Finished writing peak counts to file" << std::endl;
  }
}

std::string SpectrumFiles::getDirectory(const std::string& filepath) {
  unsigned found = filepath.find_last_of("/\\");
  return filepath.substr(0,found);
}

std::string SpectrumFiles::getFilename(const std::string& filepath) {
  unsigned found = filepath.find_last_of("/\\");
  return filepath.substr(found+1);
}

std::string SpectrumFiles::getOutputFile(const std::string& filepath, 
    const std::string& outputFolder, const std::string& newExtension) {
  std::string filename = getFilename(filepath);
  unsigned found = filename.find_last_of(".");
  return outputFolder + "/" + filename.substr(0,found) + newExtension;
}

void SpectrumFiles::createDirectory(const boost::filesystem::path& dirPath) {
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

/* MT: scanId parameter is only used in subclassed version of this class which uses feature lists */
void SpectrumFiles::getMassChargeCandidates(pwiz::msdata::SpectrumPtr s, 
    std::vector<MassChargeCandidate>& mccs, ScanId scanId) {
  SpectrumHandler::getMassChargeCandidates(s, mccs, chargeUncertainty_);
}

bool SpectrumFiles::limitsUnitTest() {
  std::string precMzFN = "/home/matthew/mergespec/data/unit_testing/Linfeng.prec_Mzs.txt";
  
  std::vector<double> precMzsAccumulated;
  std::vector<double> limits;
  
  SpectrumFiles spectrumFiles("", "");
  
  TsvInterface::read<double>(precMzFN, precMzsAccumulated);
  spectrumFiles.getPrecMzLimits(precMzsAccumulated, limits, 20.0, false);
  /*
  BOOST_FOREACH (double l, limits) {
    std::cerr << l << std::endl;
  }
  */
  int bin = -1;
  if ((bin = spectrumFiles.getPrecMzBin(precMzsAccumulated[0], limits)) != 0) {
    std::cerr << bin << " " << 0 << std::endl;
    return false;
  }
  
  if ((bin = spectrumFiles.getPrecMzBin(precMzsAccumulated[0]+1.0, limits)) != 0) {
    std::cerr << bin << " " << 0 << std::endl;
    return false;
  }
  
  if ((bin = spectrumFiles.getPrecMzBin(923.4f, limits)) != 1) {
    std::cerr << bin << " " << 0 << std::endl;
    return false;
  }
  
  if ((bin = spectrumFiles.getPrecMzBin(precMzsAccumulated[precMzsAccumulated.size()-1]-1.0, limits)) != 26) {
    std::cerr << bin << " " << 0 << std::endl;
    return false;
  }
  
  if ((bin = spectrumFiles.getPrecMzBin(precMzsAccumulated[precMzsAccumulated.size()-1], limits)) != 26) {
    std::cerr << bin << " " << 0 << std::endl;
    return false;
  }
  return true;
}

} /* namespace maracluster */
