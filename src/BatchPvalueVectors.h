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
 
#ifndef BATCH_PVALUE_VECTORS_H
#define BATCH_PVALUE_VECTORS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include "BatchGlobals.h"
#include "BatchStatement.h"
#include "BatchPvalues.h"
#include "BatchSpectrum.h"
#include "SpectrumHandler.h"
#include "SparseClustering.h"
#include "SpectrumFileList.h"
#include "PvalueCalculator.h"
#include "PvalueTriplet.h"
#include "PvalueFilterAndSort.h"
#include "PeakCounts.h"
#include "BatchPvalueVector.h"

#ifdef FINGERPRINT_FILTER
  #include "BinaryFingerprintMethods.h"
#endif

struct PvalueVectorsDbRow {
  double precMass, precMz;
  int charge;
  ScanId scannr;
  std::vector<unsigned int> peakBins;
  double retentionTime;
  int queryCharge;
  PvalueCalculator pvalCalc;
  
  inline bool operator<(const PvalueVectorsDbRow& other) const {
    return precMz < other.precMz || (precMz == other.precMz && scannr < other.scannr);
  }
};

class BatchPvalueVectors {
 public:
  BatchPvalueVectors(const std::string& pvaluesFN, double precursorTolerance, 
    bool precursorToleranceDa, double dbPvalThreshold) : 
      pvalues_(pvaluesFN), precursorTolerance_(precursorTolerance), 
      precursorToleranceDa_(precursorToleranceDa), dbPvalThreshold_(dbPvalThreshold) {}
  
  void insertMassChargeCandidate(
      MassChargeCandidate& mcc, BatchSpectrum& spec);
      
  void batchInsert(PeakCounts& peakCounts,
                   bool forceInsert);
  
  void sortPvalueVectors();
  void writePvalueVectors(const std::string& pvalueVectorsBaseFN);
  
  inline void clearPvalueVectors() { pvalVecBatch_.clear(); }
  void parsePvalueVectorFile(const std::string& pvalVecInFileFN);
  void parseBatchOverlapFile(const std::string& overlapBatchFileFN,
      std::vector< std::pair<std::string, std::string> >& overlapFNs);
  
  void processOverlapFiles(std::vector< std::pair<std::string, std::string> >& overlapFNs);
  
  void batchCalculatePvalues();
  void readFingerprints(
    std::vector<std::vector<unsigned short> >& mol_features, 
    std::vector<ScanId>& mol_identifiers, 
    std::vector<float>& prec_masses);
  void batchCalculatePvaluesJaccardFilter();
  void batchCalculatePvaluesLibrarySearch(
    std::vector<BatchSpectrum>& querySpectra);
  
  void batchCalculatePvaluesOverlap(
      std::vector<PvalueVectorsDbRow>& pvalVecCollectionTail,
      std::vector<PvalueVectorsDbRow>& pvalVecCollectionHead);
  
  static inline double getLowerBound(double precMass, double precursorTolerance, 
      bool precursorToleranceDa) {
    if (precursorToleranceDa) {
      return precMass - precursorTolerance;
    } else {
      return precMass*(1 - precursorTolerance*1e-6);
    }
  }
  static inline double getUpperBound(double precMass, double precursorTolerance, 
      bool precursorToleranceDa) {
    if (precursorToleranceDa) {
      return precMass + precursorTolerance;
    } else {
      return precMass*(1 + precursorTolerance*1e-6);
    }
  }
  static void readPvalueVectorsFile(const std::string& pvalueVectorsFN,
      std::vector<PvalueVectorsDbRow>& pvalVecCollection);
 protected:
  BatchPvalues pvalues_;
  double precursorTolerance_;
  bool precursorToleranceDa_;
  double dbPvalThreshold_;
  std::vector<PvalueVectorsDbRow> pvalVecBatch_, pvalVecCollection_;
  
  void initPvalCalc(PvalueCalculator& pvalCalc, 
                           PvalueVectorsDbRow& pvecRow, 
                           PeakCounts& peakCounts, 
                           const int numQueryPeaks, const bool polyfit);                
  
  void initPvecRow(const MassChargeCandidate& mcc, 
                          const BatchSpectrum& spec,
                          PvalueVectorsDbRow& pvecRow);
  
  void calculatePvalueVector(PvalueVectorsDbRow& pvecRow,
      PeakCounts& peakCounts);
  
  void insert(PvalueVectorsDbRow& pvecRow, 
              std::vector<BatchPvalueVector>& pvecList);
  void calculatePvalues(PvalueVectorsDbRow& pvecRow, 
                        PvalueVectorsDbRow& queryPvecRow,
                        std::vector<PvalueTriplet>& pvalBuffer);
  void calculatePvalue(PvalueVectorsDbRow& pvecRow, 
                       BatchSpectrum& querySpectrum,
                       std::vector<PvalueTriplet>& pvalBuffer);
  
  inline double getLowerBound(double mass) { 
    return getLowerBound(mass, precursorTolerance_, precursorToleranceDa_);
  }
  inline double getUpperBound(double mass) { 
    return getUpperBound(mass, precursorTolerance_, precursorToleranceDa_);
  }
  
  inline static std::string getPvalCalcKey(const std::string& uuidString, 
                                           const int charge) {
    return uuidString + "_" + boost::lexical_cast<std::string>(charge);
  }
  
  inline static bool isPvecMatch(const PvalueVectorsDbRow & pvecRow, 
                                 const PvalueVectorsDbRow & queryPvecRow) {
    return (queryPvecRow.charge == pvecRow.queryCharge) && (queryPvecRow.queryCharge == pvecRow.charge);
  }
};

#endif // BATCH_PVALUE_VECTORS_H
