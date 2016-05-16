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
 
#include "BatchSpectra.h"

void BatchSpectra::convertToBatchSpectra(SpectrumFileList& fileList) {
  std::vector<std::string> spectrumFNs = fileList.getFilePaths();
#pragma omp parallel for schedule(dynamic, 1)                
  for (int fileIdx = 0; fileIdx < spectrumFNs.size(); ++fileIdx) {
    std::string spectrumFN = spectrumFNs[fileIdx];
    convertToBatchSpectra(spectrumFN, fileList);
  }
}

void BatchSpectra::convertToBatchSpectra(std::string& spectrumFN, 
    SpectrumFileList& fileList) {
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Reading in spectra from " << spectrumFN << std::endl;
  }
  
  BatchSpectrumFiles specFiles;
  std::vector<BatchSpectrum> localSpectra;
  specFiles.getBatchSpectra(spectrumFN, fileList, localSpectra);
  
#pragma omp critical(append_spectra)
  {
    spectra_.insert(spectra_.end(), localSpectra.begin(), localSpectra.end());
  }
  
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Read " << spectra_.size() << " mass charge states." << std::endl;
  }
}

void BatchSpectra::readBatchSpectra(std::string& batchSpectraFN) {
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Reading in spectra from " << batchSpectraFN << std::endl;
  }

  if ( !boost::filesystem::exists( batchSpectraFN ) ) {
    std::cerr << "Ignoring missing file " << batchSpectraFN << std::endl;
    return;
  } 
  
  BinaryInterface::read<BatchSpectrum>(batchSpectraFN, spectra_);
  
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Read " << spectra_.size() << " mass charge states." << std::endl;
  }
}

void BatchSpectra::sortSpectraByPrecMz() {
  std::sort(spectra_.begin(), spectra_.end(), lessPrecMz);
}

void BatchSpectra::calculatePvalueVectors(SpectrumFileList& fileList, 
    PeakCounts& peakCounts) {
  if (BatchGlobals::VERB > 2) {
    std::cerr << "Inserting spectra into database" << std::endl;
  }
  
  sortSpectraByPrecMz();
  
  size_t numSpectra = spectra_.size();
  //size_t numSpectra = 20000;
  
  time_t startTime;
  time(&startTime);
  clock_t startClock = clock();

  for (size_t i = 0; i < numSpectra; ++i) {    
    if (BatchGlobals::VERB > 4) {
      std::cerr << "Global scannr " << spectra_[i].scannr << std::endl;
    }
    
    MassChargeCandidate mcc(spectra_[i].charge, spectra_[i].precMz, spectra_[i].precMass);
    pvecs_.insertMassChargeCandidate(mcc, spectra_[i]);
    
    bool forceInsert = false;
    pvecs_.batchInsert(peakCounts, forceInsert);
    
    if ((i % 50000 == 0 && BatchGlobals::VERB > 2) || BatchGlobals::VERB > 3) {
      std::cerr << "Successfully inserted spectrum " << i + 1 << "/" << 
          numSpectra << " (" << (i+1)*100/numSpectra << "%)" << std::endl;
      BatchGlobals::reportProgress(startTime, startClock, i, numSpectra);
    }
  }
  spectra_.clear();
  
  bool forceInsert = true;
  pvecs_.batchInsert(peakCounts, forceInsert);  
  pvecs_.sortPvalueVectors();
  
  if (BatchGlobals::VERB > 2) {
    std::cerr << "Successfully inserted spectra into database" << std::endl;
  }
}

void BatchSpectra::calculatePvalueVectors(PeakCounts& peakCounts) {
  SpectrumFileList fileList;
  calculatePvalueVectors(fileList, peakCounts);
}

void BatchSpectra::writePvalueVectors(const std::string& pvalueVectorsBaseFN,
                                      bool writeAll) {  
  pvecs_.writePvalueVectors(pvalueVectorsBaseFN, writeAll);
}

void BatchSpectra::calculatePvalues() {
#ifdef FINGERPRINT_FILTER
  pvecs_.batchCalculatePvaluesJaccardFilter();
#else
  pvecs_.batchCalculatePvalues();
#endif
}

void BatchSpectra::calculateAndClusterPvalues(const std::string& pvalueTreeFN,
    const std::string& scanInfoFN) {
  pvecs_.batchCalculateAndClusterPvalues(pvalueTreeFN, scanInfoFN);
}

void BatchSpectra::librarySearch(BatchSpectra& querySpectra) {
  querySpectra.sortSpectraByPrecMz();
  pvecs_.batchCalculatePvaluesLibrarySearch(querySpectra.spectra_);
}

#ifdef FINGERPRINT_FILTER
bool BatchSpectra::readFingerprints(std::string& input_file, 
    std::vector<std::vector<unsigned short> >& mol_features, 
    std::vector<ScanId>& mol_identifiers, 
    std::vector<float>& prec_masses) {  
  readBatchSpectra(input_file);
  readFingerprints(mol_features, mol_identifiers, prec_masses);
  return 1;
}

bool BatchSpectra::readFingerprints(
    std::vector<std::vector<unsigned short> >& mol_features, 
    std::vector<ScanId>& mol_identifiers, 
    std::vector<float>& prec_masses) {
  size_t numSpectra = spectra_.size();
  //size_t numSpectra = 100000;
  
  sortSpectraByPrecMz();
  
  unsigned int mol_count = 0;
  for (size_t i = 0; i < numSpectra; ++i) {
    BatchSpectrum s = spectra_[i];
    
    unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(s.precMass);
    
    std::vector<unsigned short> features;
    for (unsigned int j = 0; j < numScoringPeaks; ++j) {
      if (s.fragBins[j] != 0) {
        features.push_back(s.fragBins[j]);
      } else {
        break;
      }
    }
    if (features.size() > 0) {
      reverse(features.begin(), features.end());
      mol_identifiers.push_back(s.scannr);
      mol_features.push_back(features);
      prec_masses.push_back(s.precMass);
      ++mol_count;
    }
  }
  std::cerr << "Molecules read: " << mol_count << endl;
  return 1;
}
#endif
