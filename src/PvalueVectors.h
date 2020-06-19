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
 
#ifndef MARACLUSTER_PVALUEVECTORS_H_
#define MARACLUSTER_PVALUEVECTORS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include "Globals.h"
#include "PvalueVector.h"
#include "Pvalues.h"
#include "Spectrum.h"
#include "SpectrumFiles.h"

#include "SpectrumHandler.h"
#include "SpectrumFileList.h"

#include "PvalueCalculator.h"
#include "PvalueTriplet.h"
#include "PvalueFilterAndSort.h"

#include "PeakCounts.h"

#include "SparsePoisonedClustering.h"

#ifdef FINGERPRINT_FILTER
  #include "BinaryFingerprintMethods.h"
#endif

namespace maracluster {

struct PvalueVectorsDbRow {
  void deepCopy(const PvalueVectorsDbRow& tmp) {
    precMz = tmp.precMz;
    retentionTime = tmp.retentionTime;
    charge = tmp.charge;
    queryCharge = tmp.queryCharge;
    scannr = tmp.scannr;
    
    std::vector<unsigned int> peakBins = tmp.pvalCalc.getPeakBins();
    std::vector<unsigned int> peakScores = tmp.pvalCalc.getPeakScores();
    std::vector<double> polyfit = tmp.pvalCalc.getPolyfit();
    pvalCalc.initPolyfit(peakBins, peakScores, polyfit);
  }
  float precMz, retentionTime;
  int charge, queryCharge;
  ScanId scannr;
  PvalueCalculator pvalCalc;
  
  inline bool operator<(const PvalueVectorsDbRow& other) const {
    return precMz < other.precMz || (precMz == other.precMz && scannr < other.scannr);
  }
};

struct ClusterJob {
  size_t startBatch, endBatch;
  size_t endIdx;
  double lowerPrecMz, upperPrecMz;
  bool finished;
  std::vector<PvalueTriplet> poisonedPvals;
};

class PvalueVectors {
 public:
  PvalueVectors(const std::string& pvaluesFN, double precursorTolerance, 
    bool precursorToleranceDa, double dbPvalThreshold) : 
      pvalues_(pvaluesFN), precursorTolerance_(precursorTolerance), 
      precursorToleranceDa_(precursorToleranceDa), dbPvalThreshold_(dbPvalThreshold) {}
  
  void calculatePvalueVectors(std::vector<Spectrum>& spectra, 
      PeakCounts& peakCounts);
  
  void insertMassChargeCandidate(
      MassChargeCandidate& mcc, Spectrum& spec);
      
  void batchInsert(PeakCounts& peakCounts,
                   bool forceInsert);
  
  void sortPvalueVectors();
  void writePvalueVectors(const std::string& pvalueVectorsBaseFN, 
                          bool writeAll);
  void reloadPvalueVectors();
  inline std::vector<PvalueVectorsDbRow>& getPvalueVectors() { return pvalVecCollection_; }
  
  inline void clearPvalueVectors() { pvalVecBatch_.clear(); }
  void parsePvalueVectorFile(const std::string& pvalVecInFileFN);
  void parseBatchOverlapFile(const std::string& overlapBatchFileFN,
      std::vector< std::pair<std::string, std::string> >& overlapFNs);
  
  void processOverlapFiles(std::vector< std::pair<std::string, std::string> >& overlapFNs);
  
  void batchCalculatePvalues();
  void batchCalculateAndClusterPvalues(const std::string& resultTreeFN, 
                                       const std::string& scanInfoFN);
  void readFingerprints(
    std::vector<std::vector<unsigned short> >& mol_features, 
    std::vector<ScanId>& mol_identifiers, 
    std::vector<float>& prec_masses);
  void batchCalculatePvaluesJaccardFilter();
  void batchCalculatePvaluesLibrarySearch(
    std::vector<Spectrum>& querySpectra);
  
  void batchCalculatePvaluesOverlap(const std::string& tailFile, 
                                    const std::string& headFile);
  
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
  Pvalues pvalues_;
  double precursorTolerance_;
  bool precursorToleranceDa_;
  double dbPvalThreshold_;
  std::vector<PvalueVectorsDbRow> pvalVecBatch_, pvalVecCollection_;
  
  void initPvalCalc(PvalueCalculator& pvalCalc, 
                           PvalueVectorsDbRow& pvecRow, 
                           PeakCounts& peakCounts, 
                           const int numQueryPeaks);                
  
  void initPvecRow(const MassChargeCandidate& mcc, 
                          const Spectrum& spec,
                          PvalueVectorsDbRow& pvecRow);
  
  void calculatePvalueVector(PvalueVectorsDbRow& pvecRow,
      PeakCounts& peakCounts);
  
  void insert(PvalueVectorsDbRow& pvecRow, 
              std::vector<PvalueVector>& pvecList);
  void calculatePvalues(PvalueVectorsDbRow& pvecRow, 
                        PvalueVectorsDbRow& queryPvecRow,
                        std::vector<PvalueTriplet>& pvalBuffer);
  void calculatePvalue(PvalueVectorsDbRow& pvecRow, 
                       Spectrum& querySpectrum,
                       std::vector<PvalueTriplet>& pvalBuffer);
  
  double calculateCosineDistance(std::vector<unsigned int>& peakBins,
    std::vector<unsigned int>& queryPeakBins);
  
  void attemptClustering(size_t& newStartBatch, size_t& newPoisonedStartBatch,
    size_t pvecBatchSize, size_t numPvecBatches, size_t minPvalsForClustering,
    const std::vector<bool>& finishedPvalCalc,
    std::vector< std::vector<PvalueTriplet> >& pvalBuffers, 
    std::deque<ClusterJob>& clusterJobs, ClusterJob& poisonedClusterJob,
    std::map<ScanId, std::pair<float, float> >& precMzLimits,
    const std::string& resultTreeFN, time_t& startTime, clock_t& startClock);
    
  bool createClusterJob(size_t& newStartBatch, 
    size_t pvecBatchSize, size_t numPvecBatches, size_t minPvalsForClustering,
    const std::vector<bool>& finishedPvalCalc,
    const std::vector< std::vector<PvalueTriplet> >& pvalBuffers, 
    std::deque<ClusterJob>& clusterJobs, size_t& clusterJobIdx);
  bool createPoisonedClusterJob(size_t& newPoisonedStartBatch, 
    const size_t numPvecBatches, const size_t minPvalsForClustering, 
    std::deque<ClusterJob>& clusterJobs, ClusterJob& poisonedClusterJob);
    
  void runClusterJob(ClusterJob& clusterJob,
    std::vector< std::vector<PvalueTriplet> >& pvalBuffers,
    std::map<ScanId, std::pair<float, float> >& precMzLimits,
    const std::string& resultTreeFN,
    time_t& startTime, clock_t& startClock);
  void runPoisonedClusterJob(ClusterJob& clusterJob,
    std::deque<ClusterJob>& clusterJobs,
    std::map<ScanId, std::pair<float, float> >& precMzLimits,
    const std::string& resultTreeFN,
    const size_t numPvecBatches);
    
  void clusterPvals(std::vector<PvalueTriplet>& pvalBuffer,
    std::vector<PvalueTriplet>& pvalPoisonedBuffer,
    std::map<ScanId, std::pair<float, float> >& precMzLimits, 
    float lowerPrecMz, float upperPrecMz, const std::string& resultTreeFN);
    
  void getPrecMzLimits(
    std::map<ScanId, std::pair<float, float> >& precMzLimits);
  
  void markPoisoned(SparsePoisonedClustering& matrix, 
    std::vector<PvalueTriplet>& pvalBuffer, 
    std::map<ScanId, std::pair<float, float> >& precMzLimits, 
    float lowerPrecMz, float upperPrecMz);
  bool isPoisoned(const std::pair<float, float>& scanPrecMzLimits, 
    float lowerPrecMz, float upperPrecMz);
  bool isSafeToWrite(const std::pair<float, float>& scanPrecMzLimits, 
    float upperPrecMz);
  bool isInHead(const std::pair<float, float>& scanPrecMzLimits);
  
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

} /* namespace maracluster */

#endif /* MARACLUSTER_PVALUEVECTORS_H_ */
