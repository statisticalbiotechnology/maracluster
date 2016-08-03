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
 
#include "BatchSpectrumFiles.h"

using pwiz::msdata::SpectrumListPtr;
using pwiz::msdata::MSDataFile;
using pwiz::msdata::SpectrumPtr;

void BatchSpectrumFiles::splitByPrecursorMz(
    SpectrumFileList& fileList, std::vector<std::string>& datFNs,
    const std::string& peakCountFN, const std::string& scanInfoFN,
    double precursorTolerance, bool precursorToleranceDa) {
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Splitting spectra by precursor Mz" << std::endl;
  }
  
  std::vector<double> precMzsAccumulated;
  getPeakCountsAndPrecursorMzs(fileList, precMzsAccumulated, peakCountFN);
  
  //writePrecMzs(precMzsAccumulated);
  //readPrecMzs(peakCountFN_, precMzsAccumulated);
  
  std::vector<double> limits;
  getPrecMzLimits(precMzsAccumulated, limits, precursorTolerance, 
                    precursorToleranceDa);
  getDatFNs(limits, datFNs);
  writeSplittedPrecursorMzFiles(fileList, limits, datFNs, scanInfoFN);
}

void BatchSpectrumFiles::splitByPrecursorMz(SpectrumFileList& fileList,
    const std::string& datFNFile, const std::string& peakCountFN,
    const std::string& scanInfoFN, double precursorTolerance, 
    bool precursorToleranceDa) {  
  std::vector<std::string> datFNs;
  splitByPrecursorMz(fileList, datFNs, peakCountFN, scanInfoFN, 
                       precursorTolerance, precursorToleranceDa);
  
  writeDatFNsToFile(datFNs, datFNFile);
}

void BatchSpectrumFiles::getPeakCountsAndPrecursorMzs(
    SpectrumFileList& fileList, 
    std::vector<double>& precMzsAccumulated,
    const std::string& peakCountFN) {
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Accumulating peak counts and precursor Mzs" << std::endl;
  }
  
  PeakCounts peakCountsAccumulated;
  
  std::vector<std::string> spectrumFNs = fileList.getFilePaths();
#pragma omp parallel for schedule(dynamic, 1)  
  for (int fileIdx = 0; fileIdx < static_cast<int>(spectrumFNs.size()); ++fileIdx) {
    std::string spectrumFN = spectrumFNs[fileIdx];
    if (BatchGlobals::VERB > 1) {
      std::cerr << "  Processing " << spectrumFN << 
          " (" << (fileIdx+1)*100/spectrumFNs.size() << "%)." << std::endl;
    }
    
    SpectrumListPtr specList;    
  #pragma omp critical (create_msdata)
    {  
      MSReaderList readerList;
      MSDataFile msd(spectrumFN, &readerList);
      specList = msd.run.spectrumListPtr;
    }
    PeakCounts peakCounts;
    std::vector<double> precMzs;
    
    size_t numSpectra = specList->size();
    //size_t numSpectra = 2;
    for (size_t i = 0; i < numSpectra; ++i) {
      SpectrumPtr s = specList->spectrum(i, true);
      if (!SpectrumHandler::isMs2Scan(s)) continue;
      
      std::vector<MZIntensityPair> mziPairs;
      SpectrumHandler::getMZIntensityPairs(s, mziPairs); 
      
      std::vector<MassChargeCandidate> mccs;
      SpectrumHandler::getMassChargeCandidates(s, mccs, chargeUncertainty_); // returns mccs sorted by charge
      unsigned int lastCharge = 0;
      BOOST_FOREACH (MassChargeCandidate& mcc, mccs) {
        precMzs.push_back(mcc.precMz);
        unsigned int charge = (std::min)(mcc.charge, peakCounts.getMaxCharge());
        if (charge != lastCharge) {
          // in the last bin we do not truncate the spectrum
          if (charge == peakCounts.getMaxCharge()) charge = 100u;
          unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(mcc.mass);
          peakCounts.addSpectrum(mziPairs, mcc.precMz, charge, mcc.mass, numScoringPeaks);
          lastCharge = charge;
        }
      }
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

void BatchSpectrumFiles::writeSplittedPrecursorMzFiles(
    SpectrumFileList& fileList, 
    std::vector<double>& limits,
    std::vector<std::string>& datFNs,
    const std::string& scanInfoFN) {
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Dividing spectra in " << limits.size() << 
                 " bins of ~2 CPU hours each." << std::endl;
  }
  
  std::vector<std::string> spectrumFNs = fileList.getFilePaths();
#pragma omp parallel for schedule(dynamic, 1)                
  for (int fileIdx = 0; fileIdx < static_cast<int>(spectrumFNs.size()); ++fileIdx) {
    std::string spectrumFN = spectrumFNs[fileIdx];
    if (BatchGlobals::VERB > 1) {
      std::cerr << "  Processing " << spectrumFN << 
          " (" << (fileIdx+1)*100/spectrumFNs.size() << "%)." << std::endl;
    }
    
    std::vector<BatchSpectrum> localSpectra;
    std::vector<ScanInfo> scanInfos;
    getBatchSpectra(spectrumFN, fileList, localSpectra, scanInfos);
    
    std::vector< std::vector<BatchSpectrum> > batchSpectra(limits.size());
    BOOST_FOREACH (BatchSpectrum& bs, localSpectra) {
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
    }
  }
}


void BatchSpectrumFiles::getBatchSpectra(
    const std::string& spectrumFN, SpectrumFileList& fileList,
    std::vector<BatchSpectrum>& localSpectra, 
    std::vector<ScanInfo>& scanInfos) {
  if ( !boost::filesystem::exists( spectrumFN ) ) {
    std::cerr << "Ignoring missing file " << spectrumFN << std::endl;
    return;
  }
   
  SpectrumListPtr specList;
#pragma omp critical (create_msdata)
  {
    MSReaderList readerList;
    MSDataFile msd(spectrumFN, &readerList);
    specList = msd.run.spectrumListPtr;
  }
  
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
    
    std::vector<MassChargeCandidate> mccs;
    SpectrumHandler::getMassChargeCandidates(s, mccs, chargeUncertainty_);
    
    BOOST_FOREACH (MassChargeCandidate& mcc, mccs) {
      for (int isotopeTolerance = 0; isotopeTolerance <= 0; ++isotopeTolerance) {
        double mass = mcc.mass + isotopeTolerance;
        std::vector<unsigned int> peakBins;
        unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(mass);
        BinSpectra::binBinaryTruncated(mziPairs, peakBins, 
          numScoringPeaks, mcc.mass);
        
        if (peakBins.size() >= PvalueCalculator::getMinScoringPeaks(mass)) {
          BatchSpectrum bs;
          peakBins.resize(BATCH_SPECTRUM_NUM_STORED_PEAKS, 0u);
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

void BatchSpectrumFiles::appendBatchSpectra(
    std::vector< std::vector<BatchSpectrum> >& batchSpectra,
    std::vector<std::string>& datFNs) {
  bool append = true;
  for (size_t i = 0; i < datFNs.size(); ++i) {
    BinaryInterface::write<BatchSpectrum>(batchSpectra[i], datFNs[i], append);
  }
}

void BatchSpectrumFiles::getDatFNs(std::vector<double>& limits, 
    std::vector<std::string>& datFNs) {
  int lastLimit = -1;
  int counter = 0;
  BOOST_FOREACH (double limit, limits) {
    int intLimit = static_cast<int>(limit);
    std::string filepath = precMzFileFolder_ + "/" + 
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

void BatchSpectrumFiles::readDatFNsFromFile(const std::string& datFNFile,
    std::vector<std::string>& datFNs) {
  std::ifstream infile(datFNFile.c_str(), std::ios_base::in);
  if (infile.is_open()) {
    std::string datFN;
    while (getline(infile, datFN)) {
      datFNs.push_back(datFN);
    }
  } else {
    std::cerr << "Could not read list of dat files" << std::endl;
  }
}

void BatchSpectrumFiles::writeDatFNsToFile(std::vector<std::string>& datFNs,
    const std::string& datFNFile) {
  std::ofstream outfile(datFNFile.c_str(), std::ios_base::out);
  if (outfile.is_open()) {
    BOOST_FOREACH (std::string& datFN, datFNs) {
      outfile << datFN << "\n";
    }
  } else {
    std::cerr << "Could not write list of dat files" << std::endl;
  }
}

void BatchSpectrumFiles::writePrecMzs(const std::vector<double>& precMzs) {
  BOOST_FOREACH (const double precMz, precMzs) {
    std::cout << precMz << std::endl;
  }
}

void BatchSpectrumFiles::readPrecMzs(const std::string& precMzFN,
    std::vector<double>& precMzs) {
  std::ifstream precMzStream;
  precMzStream.open(precMzFN.c_str(), std::ios::in | std::ios::binary);
  if (!precMzStream.is_open()) {
    std::cerr << "Could not open file " << precMzFN << std::endl;
  } else {
    std::string line;
    while (getline(precMzStream, line)) {
      precMzs.push_back(boost::lexical_cast<double>(line));
    }
  }
}

void BatchSpectrumFiles::readPrecMzLimits(const std::string& scanInfoFN,
    std::map<ScanId, std::pair<float, float> >& precMzLimits) {
  if (BatchGlobals::VERB > 2) {
    std::cerr << "Reading precursor m/z limits." << std::endl;
  }
  std::vector<ScanInfo> scanInfos;
  BinaryInterface::read<ScanInfo>(scanInfoFN, scanInfos);
  
  BOOST_FOREACH (const ScanInfo& si, scanInfos) {
    precMzLimits[si.scanId] = std::make_pair(si.minPrecMz, si.maxPrecMz);
  }
}

void BatchSpectrumFiles::getPrecMzLimits(std::vector<double>& precMzs, 
    std::vector<double>& limits, double precursorTolerance, 
    bool precursorToleranceDa) {
  float pvecCost = 0.7f; // computation time for 1 p-value vector in ms on 1 core
  float pvalCost = 0.001f; // computation time for 1 p-value pair in ms on 1 core
  float maxCost = 120.0f*60.0f*1000.0f; // max computation time = 120 CPU minutes
  
  unsigned long long numComparisons = 0uL;
  if (precMzs.size() > 0) {
    size_t lowerBoundIdx = 0;
    limits.push_back(precMzs[0]);
    float curCost = 0.0f;
    for (size_t idx = 0; idx < precMzs.size(); ++idx) {
      curCost += pvecCost;
      double lowerBound = BatchPvalueVectors::getLowerBound(precMzs[idx], 
          precursorTolerance, precursorToleranceDa);
      while (precMzs[lowerBoundIdx] < lowerBound) {
        ++lowerBoundIdx;
      }
      numComparisons += (idx - lowerBoundIdx);
      curCost += (idx - lowerBoundIdx)*pvalCost;
      if (curCost > maxCost && lowerBound > limits.back()) {
        curCost = 0.0f;
        limits.push_back(precMzs[idx]);
      }
    } 
  }
  
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Estimated number of pair comparisons: " 
              << numComparisons << std::endl;
  }
}

int BatchSpectrumFiles::getPrecMzBin(double precMz, std::vector<double>& limits) {
  int bin = std::upper_bound(limits.begin(), limits.end(), precMz) - limits.begin() - 1;
  return std::min(std::max(bin, 0), static_cast<int>(limits.size()) - 1);
}

void BatchSpectrumFiles::writePeakCounts(PeakCounts& peakCountsAccumulated, 
                                         const std::string& peakCountFN) {
  if (BatchGlobals::VERB > 2) {
    std::cerr << "Writing peak counts to file" << std::endl;
  }
  std::string serializedPeakCounts;
  PeakCounts::serializePeakCounts(peakCountsAccumulated, serializedPeakCounts);
  
  std::ofstream peakCountStream(peakCountFN.c_str(), std::ios_base::binary | std::ios_base::out);
  if (peakCountStream.is_open()) {
    peakCountStream << serializedPeakCounts;
  }
  
  if (BatchGlobals::VERB > 2) {
    std::cerr << "Finished writing peak counts to file" << std::endl;
  }
}

std::string BatchSpectrumFiles::getDirectory(const std::string& filepath) {
  unsigned found = filepath.find_last_of("/\\");
  return filepath.substr(0,found);
}

std::string BatchSpectrumFiles::getFilename(const std::string& filepath) {
  unsigned found = filepath.find_last_of("/\\");
  return filepath.substr(found+1);
}

bool BatchSpectrumFiles::limitsUnitTest() {
  std::string precMzFN = "/home/matthew/mergespec/data/unit_testing/Linfeng.prec_Mzs.txt";
  
  std::vector<double> precMzsAccumulated;
  std::vector<double> limits;
  
  BatchSpectrumFiles spectrumFiles("");
  
  spectrumFiles.readPrecMzs(precMzFN, precMzsAccumulated);
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
