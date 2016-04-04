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
using pwiz::msdata::SpectrumListSimple;
using pwiz::msdata::SpectrumListSimplePtr;
using pwiz::msdata::MSData;
using pwiz::msdata::MSDataMerger;
using pwiz::msdata::MSDataPtr;
using pwiz::msdata::MSDataFile;
using pwiz::msdata::SpectrumPtr;
using pwiz::msdata::Spectrum;
using pwiz::msdata::SelectedIon;

void BatchSpectrumFiles::splitByPrecursorMass(
    SpectrumFileList& fileList, std::vector<std::string>& datFNs,
    const std::string& peakCountFN, const std::string& scanNrsFN,
    double precursorTolerance, bool precursorToleranceDa) {
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Splitting spectra by precursor mass" << std::endl;
  }
  
  std::vector<double> precMassesAccumulated;
  getPeakCountsAndPrecursorMasses(fileList, precMassesAccumulated, peakCountFN);
  
  //writePrecMasses(precMassesAccumulated);
  //readPrecMasses(peakCountFN_, precMassesAccumulated);
  
  std::vector<double> limits;
  getPrecMassLimits(precMassesAccumulated, limits, precursorTolerance, 
                    precursorToleranceDa);
  getDatFNs(limits, datFNs);
  writeSplittedPrecursorMassFiles(fileList, limits, datFNs, scanNrsFN);
}

void BatchSpectrumFiles::splitByPrecursorMass(SpectrumFileList& fileList,
    const std::string& datFNFile, const std::string& peakCountFN,
    const std::string& scanNrsFN, double precursorTolerance, 
    bool precursorToleranceDa) {  
  std::vector<std::string> datFNs;
  splitByPrecursorMass(fileList, datFNs, peakCountFN, scanNrsFN, 
                       precursorTolerance, precursorToleranceDa);
  
  writeDatFNsToFile(datFNs, datFNFile);
}

void BatchSpectrumFiles::getPeakCountsAndPrecursorMasses(
    SpectrumFileList& fileList, 
    std::vector<double>& precMassesAccumulated,
    const std::string& peakCountFN) {
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Accumulating peak counts and precursor masses" << std::endl;
  }
  
  PeakCounts peakCountsAccumulated;
  
  std::vector<std::string> spectrumFNs = fileList.getFilePaths();
#pragma omp parallel for schedule(dynamic, 1)  
  for (int fileIdx = 0; fileIdx < spectrumFNs.size(); ++fileIdx) {
    std::string spectrumFN = spectrumFNs[fileIdx];
    if (BatchGlobals::VERB > 1) {
      std::cerr << "  Processing " << spectrumFN << 
          " (" << (fileIdx+1)*100/spectrumFNs.size() << "%)." << std::endl;
    }
    
    SpectrumListPtr specList;    
  #pragma omp critical (create_msdata)
    {  
      MSDataFile msd(spectrumFN);
      specList = msd.run.spectrumListPtr;
    }
    PeakCounts peakCounts;
    std::vector<double> precMasses;
    
    size_t numSpectra = specList->size();
    //size_t numSpectra = 2;
    
    for (size_t i = 0; i < numSpectra; ++i) {
      SpectrumPtr s = specList->spectrum(i, true);
      
      std::vector<MZIntensityPair> mziPairs;
      SpectrumHandler::getMZIntensityPairs(s, mziPairs); 
      
      std::vector<MassChargeCandidate> mccs;
      SpectrumHandler::getMassChargeCandidates(s, mccs, chargeUncertainty_); // returns mccs sorted by charge
      unsigned int lastCharge = 0;
      BOOST_FOREACH (MassChargeCandidate& mcc, mccs) {
        precMasses.push_back(mcc.mass);
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
      precMassesAccumulated.insert( precMassesAccumulated.end(), precMasses.begin(), precMasses.end() );
    }
  }
  writePeakCounts(peakCountsAccumulated, peakCountFN);
  
  std::sort(precMassesAccumulated.begin(), precMassesAccumulated.end());
}

void BatchSpectrumFiles::writeScannrs(SpectrumFileList& fileList,
                                      const std::string& scanNrsFN) {
  std::vector<std::string> spectrumFNs = fileList.getFilePaths();
#pragma omp parallel for schedule(dynamic, 1)                
  for (int fileIdx = 0; fileIdx < spectrumFNs.size(); ++fileIdx) {
    std::string spectrumFN = spectrumFNs[fileIdx];
    
    SpectrumListPtr specList;    
  #pragma omp critical (create_msdata)
    {  
      MSDataFile msd(spectrumFN);
      specList = msd.run.spectrumListPtr;
    }
    size_t numSpectra = specList->size();
    //size_t numSpectra = 2;
    std::vector<ScanId> globalScanNrs(numSpectra);
    
    for(size_t i = 0; i < numSpectra; ++i) {
      SpectrumPtr s = specList->spectrum(i, false);
      
      unsigned int scannr = SpectrumHandler::getScannr(s);
      ScanId globalIdx = fileList.getScanId(spectrumFN, scannr);
      globalScanNrs[i] = globalIdx;
    }
    bool append = true;
    BinaryInterface::write<ScanId>(globalScanNrs, scanNrsFN, append);
  }
}

void BatchSpectrumFiles::writeSplittedPrecursorMassFiles(
    SpectrumFileList& fileList, 
    std::vector<double>& limits,
    std::vector<std::string>& datFNs,
    const std::string& scanNrsFN) {
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Dividing spectra in " << limits.size() << 
                 " bins of ~2 CPU hours each." << std::endl;
  }
  
  std::vector<std::string> spectrumFNs = fileList.getFilePaths();
#pragma omp parallel for schedule(dynamic, 1)                
  for (int fileIdx = 0; fileIdx < spectrumFNs.size(); ++fileIdx) {
    std::string spectrumFN = spectrumFNs[fileIdx];
    if (BatchGlobals::VERB > 1) {
      std::cerr << "  Processing " << spectrumFN << 
          " (" << (fileIdx+1)*100/spectrumFNs.size() << "%)." << std::endl;
    }
    
    SpectrumListPtr specList;    
  #pragma omp critical (create_msdata)
    {  
      MSDataFile msd(spectrumFN);
      specList = msd.run.spectrumListPtr;
    }
    
    std::vector< std::vector<BatchSpectrum> > batchSpectra(limits.size());
    size_t numSpectra = specList->size();
    //size_t numSpectra = 2;
    std::vector<ScanId> globalScanNrs(numSpectra);
    
    for (size_t i = 0; i < numSpectra; ++i) {
      SpectrumPtr s = specList->spectrum(i, true);
      
      std::vector<MZIntensityPair> mziPairs;
      SpectrumHandler::getMZIntensityPairs(s, mziPairs); 
      
      double retentionTime = SpectrumHandler::getRetentionTime(s);
      unsigned int scannr = SpectrumHandler::getScannr(s);
      ScanId globalIdx = fileList.getScanId(spectrumFN, scannr);
      
      globalScanNrs[i] = globalIdx;
      
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
            bs.precMass = mass;
            bs.retentionTime = retentionTime;
            bs.charge = mcc.charge;
            bs.scannr = globalIdx;
            
            int precBin = getPrecMassBin(mass, limits);
            batchSpectra[precBin].push_back(bs);
          }
        }
      }
    }
  #pragma omp critical (add_to_datfiles)
    {
      appendBatchSpectra(batchSpectra, datFNs);
    }
  #pragma omp critical (write_scannrs)
    {
      bool append = true;
      BinaryInterface::write<ScanId>(globalScanNrs, scanNrsFN, append);
    }
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
    std::string filepath = precMassFileFolder_ + "/" + 
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

void BatchSpectrumFiles::writePrecMasses(const std::vector<double>& precMasses) {
  BOOST_FOREACH (const double precMass, precMasses) {
    std::cout << precMass << std::endl;
  }
}

void BatchSpectrumFiles::readPrecMasses(const std::string& precMassFN,
    std::vector<double>& precMasses) {
  std::ifstream precMassStream;
  precMassStream.open(precMassFN.c_str(), std::ios::in | std::ios::binary);
  if (!precMassStream.is_open()) {
    std::cerr << "Could not open file " << precMassFN << std::endl;
  } else {
    std::string line;
    while (getline(precMassStream, line)) {
      precMasses.push_back(boost::lexical_cast<double>(line));
    }
  }
}

void BatchSpectrumFiles::getPrecMassLimits(std::vector<double>& precMasses, 
    std::vector<double>& limits, double precursorTolerance, 
    bool precursorToleranceDa) {
  float pvecCost = 0.7f; // computation time for 1 p-value vector in ms on 1 core
  float pvalCost = 0.001f; // computation time for 1 p-value pair in ms on 1 core
  float maxCost = 120.0f*60.0f*1000.0f; // max computation time = 120 CPU minutes
  float curCost = 0.0f;
  
  unsigned long long numComparisons = 0uL;
  
  if (precMasses.size() > 0) {
    size_t lowerBoundIdx = 0;
    limits.push_back(precMasses[0]);
    for (size_t idx = 0; idx < precMasses.size(); ++idx) {
      curCost += pvecCost;
      double lowerBound = BatchPvalueVectors::getLowerBound(precMasses[idx], 
          precursorTolerance, precursorToleranceDa);
      while (precMasses[lowerBoundIdx] < lowerBound) {
        ++lowerBoundIdx;
      }
      numComparisons += (idx - lowerBoundIdx);
      curCost += (idx - lowerBoundIdx)*pvalCost;
      if (curCost > maxCost) {
        curCost = 0.0f;
        limits.push_back(precMasses[idx]);
      }
    } 
  }
  
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Estimated number of pair comparisons: " 
              << numComparisons << std::endl;
  }
}

int BatchSpectrumFiles::getPrecMassBin(double precMass, std::vector<double>& limits) {
  int bin = std::upper_bound(limits.begin(), limits.end(), precMass) - limits.begin() - 1;
  return std::min(std::max(bin, 0), static_cast<int>(limits.size()) - 1);
}

void BatchSpectrumFiles::writePeakCounts(PeakCounts& peakCountsAccumulated, 
                                         const std::string& peakCountFN) {
  if (BatchGlobals::VERB > 2) {
    std::cerr << "Writing peak counts to file" << std::endl;
  }
  std::string serializedPeakCounts;
  PeakCounts::serializePeakCounts(peakCountsAccumulated, serializedPeakCounts);
  
  std::ofstream peakCountStream(peakCountFN.c_str());
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
  std::string precMassFN = "/home/matthew/mergespec/data/unit_testing/Linfeng.prec_masses.txt";
  
  std::vector<double> precMassesAccumulated;
  std::vector<double> limits;
  
  BatchSpectrumFiles spectrumFiles("");
  
  spectrumFiles.readPrecMasses(precMassFN, precMassesAccumulated);
  spectrumFiles.getPrecMassLimits(precMassesAccumulated, limits, 20.0, false);
  /*
  BOOST_FOREACH (double l, limits) {
    std::cerr << l << std::endl;
  }
  */
  int bin = -1;
  if ((bin = spectrumFiles.getPrecMassBin(precMassesAccumulated[0], limits)) != 0) {
    std::cerr << bin << " " << 0 << std::endl;
    return false;
  }
  
  if ((bin = spectrumFiles.getPrecMassBin(precMassesAccumulated[0]+1.0, limits)) != 0) {
    std::cerr << bin << " " << 0 << std::endl;
    return false;
  }
  
  if ((bin = spectrumFiles.getPrecMassBin(923.4f, limits)) != 1) {
    std::cerr << bin << " " << 0 << std::endl;
    return false;
  }
  
  if ((bin = spectrumFiles.getPrecMassBin(precMassesAccumulated[precMassesAccumulated.size()-1]-1.0, limits)) != 26) {
    std::cerr << bin << " " << 0 << std::endl;
    return false;
  }
  
  if ((bin = spectrumFiles.getPrecMassBin(precMassesAccumulated[precMassesAccumulated.size()-1], limits)) != 26) {
    std::cerr << bin << " " << 0 << std::endl;
    return false;
  }
  return true;
}
