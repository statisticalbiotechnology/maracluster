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
 
#ifndef BATCH_SPECTRA_H
#define BATCH_SPECTRA_H

#include <iostream>
#include <vector>
#include <string>

#include <boost/filesystem.hpp>

#include "pwiz/data/msdata/MSDataFile.hpp"

#include "BatchGlobals.h"
#include "BatchSpectrumFiles.h"
#include "BatchSpectrum.h"
#include "BatchPvalueVectors.h"

#include "SpectrumHandler.h"
#include "SpectrumFileList.h"
#include "PvalueCalculator.h"
#include "PeakCounts.h"

class BatchSpectra {
 public:
  BatchSpectra(const std::string& pvaluesFN) : pvecs_(pvaluesFN) {}
  
  // methods to import spectra
  void setBatchSpectra(std::vector<BatchSpectrum>& spectra) {
    spectra_ = spectra;
  }
  void convertToBatchSpectra(std::string& spectrumFN, 
    SpectrumFileList& fileList);
  void convertToBatchSpectra(SpectrumFileList& fileList);
  void readBatchSpectra(std::string& batchSpectraFN);
  
  // methods for p-value vectors
  void calculatePvalueVectors(SpectrumFileList& fileList, 
    PeakCounts& peakCounts);
  void calculatePvalueVectors(PeakCounts& peakCounts);
  void writePvalueVectors(const std::string& pvalueVectorsBaseFN);
  
  // method to calculate p-values
  void calculatePvalues();
  void librarySearch(BatchSpectra& querySpectra);
  
  // Reading a batchspectrum inputfile for fingerprint similarities
  bool readFingerprints(std::string& input_file, 
      std::vector<std::vector<unsigned short> >& mol_features, 
      std::vector<ScanId>& mol_identifiers, 
      std::vector<float>& prec_masses);
  bool readFingerprints(
      std::vector<std::vector<unsigned short> >& mol_features, 
      std::vector<ScanId>& mol_identifiers, 
      std::vector<float>& prec_masses);
  
  inline static bool lessPrecMass(const BatchSpectrum& a, 
    const BatchSpectrum& b) { return (a.precMass < b.precMass) || (a.precMass == b.precMass && a.scannr < b.scannr); }
 protected:
  BatchPvalueVectors pvecs_;
  std::vector<BatchSpectrum> spectra_;
  
  void sortSpectraByPrecMass();
};

#endif // BATCH_SPECTRA_H
