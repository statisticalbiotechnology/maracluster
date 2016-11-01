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
 
#include "BatchPvalueVectors.h"

void BatchPvalueVectors::calculatePvalueVectors(
    std::vector<BatchSpectrum>& spectra, 
    PeakCounts& peakCounts) {
  if (Globals::VERB > 2) {
    std::cerr << "Inserting spectra into database" << std::endl;
  }
  
  size_t numSpectra = spectra.size();
  //size_t numSpectra = 500000;
  
  time_t startTime;
  time(&startTime);
  clock_t startClock = clock();
  
  for (size_t i = 0; i < numSpectra; ++i) {    
    if (Globals::VERB > 4) {
      std::cerr << "Global scannr " << spectra[i].scannr << std::endl;
    }
    
    float precMass = SpectrumHandler::calcMass(spectra[i].precMz, 
                                               spectra[i].charge);
    MassChargeCandidate mcc(spectra[i].charge, spectra[i].precMz, precMass);
    insertMassChargeCandidate(mcc, spectra[i]);
    
    bool forceInsert = false;
    batchInsert(peakCounts, forceInsert);
    
    if ((i % 50000 == 0 && Globals::VERB > 2) || Globals::VERB > 3) {
      std::cerr << "Successfully inserted spectrum " << i + 1 << "/" << 
          numSpectra << " (" << (i+1)*100/numSpectra << "%)" << std::endl;
      Globals::reportProgress(startTime, startClock, i, numSpectra);
    }
  }
  
  bool forceInsert = true;
  batchInsert(peakCounts, forceInsert);  
  sortPvalueVectors();
  reloadPvalueVectors();
  
  if (Globals::VERB > 2) {
    std::cerr << "Successfully inserted spectra into database" << std::endl;
  }
}

void BatchPvalueVectors::initPvecRow(const MassChargeCandidate& mcc, 
                                    const BatchSpectrum& spec,
                                    PvalueVectorsDbRow& pvecRow) {
  pvecRow.precMz = mcc.precMz;
  pvecRow.charge = mcc.charge;
  pvecRow.scannr = spec.scannr;
  pvecRow.retentionTime = spec.retentionTime;
  
  float precMass = SpectrumHandler::calcMass(spec.precMz, spec.charge);
  unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(precMass);
  std::vector<unsigned int> peakBins;
#ifdef DOT_PRODUCT
  numScoringPeaks *= 2;
  peakBins.reserve(numScoringPeaks);
  for (unsigned int j = 0; j < numScoringPeaks; ++j) {
    if (spec.fragBins[j] != 0 || j % 2 == 1) {
      peakBins.push_back(spec.fragBins[j]);
    } else if (j % 2 == 0) {
      break;
    }
  }
#else
  peakBins.reserve(numScoringPeaks);
  for (unsigned int j = 0; j < numScoringPeaks; ++j) {
    if (spec.fragBins[j] != 0) {
      peakBins.push_back(spec.fragBins[j]);
    } else {
      break;
    }
  }
#endif
  pvecRow.pvalCalc.setPeakBins(peakBins);
}

void BatchPvalueVectors::insertMassChargeCandidate(
    MassChargeCandidate& mcc, BatchSpectrum& spec) {
  if (Globals::VERB > 4) {
    std::cerr << "Inserting mass charge candidate into database" << std::endl;
  }
  
  /*
  int minCharge = (std::max)(mcc.charge - 1, 1u);
  int maxCharge = mcc.charge + 1;
  */
  int minCharge = mcc.charge;
  int maxCharge = mcc.charge;
  for (int charge = minCharge; charge <= maxCharge; ++charge) {
    PvalueVectorsDbRow pvecRow;
    initPvecRow(mcc, spec, pvecRow);
    
    pvecRow.queryCharge = charge;
    pvalVecBatch_.push_back(pvecRow);
  }
  
  if (Globals::VERB > 4) {
    std::cerr << "Inserted mass charge candidate into database" << std::endl;
  }
}

void BatchPvalueVectors::batchInsert(PeakCounts& peakCounts, bool forceInsert) {
  if (pvalVecBatch_.size() >= 5000 || forceInsert) {
    if (Globals::VERB > 3) {
      std::cerr << "Inserting " << pvalVecBatch_.size() << " spectra into database" << std::endl;
    }
    //std::cerr << "Calculating " << pvalVecBatch_.size() << " p-value vectors" << std::endl;
  #pragma omp parallel for schedule(dynamic, 100)
    for (int i = 0; i < pvalVecBatch_.size(); ++i) {
      calculatePvalueVector(pvalVecBatch_[i], peakCounts);
    }
    pvalVecBatch_.clear();
    
    if (Globals::VERB > 3) {
      std::cerr << "Currently there are " << pvalVecCollection_.size() << " spectra in the database" << std::endl;
    }
  }
}

void BatchPvalueVectors::initPvalCalc(PvalueCalculator& pvalCalc, 
    PvalueVectorsDbRow& pvecRow, PeakCounts& peakCounts, 
    const int numQueryPeaks) {
#ifndef DOT_PRODUCT
  PeakDistribution distribution;
  peakCounts.generatePeakDistribution(pvecRow.precMz, pvecRow.queryCharge, 
                                      distribution, numQueryPeaks);
  
  pvalCalc.computePvalVectorPolyfit(distribution.getDistribution());
#endif
}

void BatchPvalueVectors::calculatePvalueVector(PvalueVectorsDbRow& pvecRow,
    PeakCounts& peakCounts) {
  if (Globals::VERB > 4) {
    std::cerr << "Inserting pvalue vector " << pvecRow.scannr << std::endl;
  }
  float precMass = SpectrumHandler::calcMass(pvecRow.precMz, pvecRow.charge);
  unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(precMass);
  initPvalCalc(pvecRow.pvalCalc, pvecRow, peakCounts, numScoringPeaks);
  
  if (pvecRow.pvalCalc.getNumScoringPeaks() >= PvalueCalculator::getMinScoringPeaks(precMass)) {    
  #pragma omp critical (store_pvec)
    {
      pvalVecCollection_.push_back(pvecRow);
    }
        
    if (Globals::VERB > 4) {
      std::cerr << "Inserted pvalue vector " << pvecRow.scannr << std::endl;
    }
  } else if (Globals::VERB > 2) {
    std::cerr << "Ignoring " << pvecRow.scannr << ". Not enough scoring peaks: " << pvecRow.pvalCalc.getNumScoringPeaks() << std::endl;
  }
}

void BatchPvalueVectors::sortPvalueVectors() {
  if (Globals::VERB > 2) {
    std::cerr << "Sorting pvalue vectors" << std::endl;
  }
  std::sort(pvalVecCollection_.begin(), pvalVecCollection_.end());
}

void BatchPvalueVectors::writePvalueVectors(
    const std::string& pvalueVectorsBaseFN, bool writeAll) {
  std::string pvalueVectorsFN = pvalueVectorsBaseFN + ".dat";
  std::string pvalueVectorsHeadFN = pvalueVectorsBaseFN + ".head.dat";
  std::string pvalueVectorsTailFN = pvalueVectorsBaseFN + ".tail.dat";
  
  if (Globals::VERB > 1) {
    std::cerr << "Writing pvalue vectors" << std::endl;
  }

  std::vector<BatchPvalueVector> headList, tailList, allList;
  double headOverlapLimit = getUpperBound(pvalVecCollection_.front().precMz);
  double tailOverlapLimit = getLowerBound(pvalVecCollection_.back().precMz);
  size_t n = pvalVecCollection_.size();
  
  for (size_t i = 0; i < n; ++i) {
    if (i % 100000 == 0 && Globals::VERB > 2) {
      std::cerr << "Writing pvalue vector " << i << "/" << n << std::endl;
    }

    if (writeAll) insert(pvalVecCollection_[i], allList);
    
    if (pvalVecCollection_[i].precMz < headOverlapLimit) {
      insert(pvalVecCollection_[i], headList);
    }
    if (pvalVecCollection_[i].precMz > tailOverlapLimit) {
      insert(pvalVecCollection_[i], tailList);
    }
  }

  bool append = false;
  BinaryInterface::write<BatchPvalueVector>(allList, pvalueVectorsFN, append);
  BinaryInterface::write<BatchPvalueVector>(headList, pvalueVectorsHeadFN, append);
  BinaryInterface::write<BatchPvalueVector>(tailList, pvalueVectorsTailFN, append);
  
  if (Globals::VERB > 1) {
    std::cerr << "Finished writing pvalue vectors" << std::endl;
  }
}

void BatchPvalueVectors::insert(PvalueVectorsDbRow& pvecRow, std::vector<BatchPvalueVector>& pvecList) {
  if (Globals::VERB > 4) {
    std::cerr << "Inserting pvalue vector into pvalue vectors " <<
                 "table asynchronously" << std::endl;
  }  
  BatchPvalueVector pvec;
  
  pvec.precMz = pvecRow.precMz;
  pvec.charge = pvecRow.charge;
  pvec.scannr = pvecRow.scannr;
  
  pvec.retentionTime = pvecRow.retentionTime;
  pvec.queryCharge = pvecRow.queryCharge;
  
  pvecRow.pvalCalc.copyPolyfit(pvec.peakBins, pvec.peakScores, pvec.polyfit);
  
  pvecList.push_back(pvec);
  
  if (Globals::VERB > 4) {
    std::cerr << "Put pvalue vector insertion into queue" << std::endl;
  }
}

void BatchPvalueVectors::reloadPvalueVectors() {
  std::vector<PvalueVectorsDbRow> pvalVecCollection;
  BOOST_FOREACH (PvalueVectorsDbRow& tmp, pvalVecCollection_) {
    pvalVecCollection.push_back(tmp);
  }
  pvalVecCollection_.swap(pvalVecCollection);
  
  if (Globals::VERB > 1) {
    std::cerr << "Reloaded " << pvalVecCollection.size() << " pvalue vectors." << std::endl;
  }
}

/*
void BatchPvalueVectors::transferPvalueVectors(BatchPvalueVectors& pvecs) {
  BOOST_FOREACH (PvalueVectorsDbRow& tmp, pvecs.getPvalueVectors()) {
    pvalVecCollection_.push_back(tmp);
  }
}
*/

void BatchPvalueVectors::readPvalueVectorsFile(const std::string& pvalueVectorsFN,
    std::vector<PvalueVectorsDbRow>& pvalVecCollection) {  
  if (Globals::VERB > 1) {
    std::cerr << "Reading in pvalue vectors from " << pvalueVectorsFN << std::endl;
  }

  if (Globals::fileIsEmpty(pvalueVectorsFN)) {
    std::cerr << "Ignoring missing/empty file " << pvalueVectorsFN << std::endl;
    return;
  }

  boost::iostreams::mapped_file mmap(pvalueVectorsFN,
            boost::iostreams::mapped_file::readonly);
  
  const char* f = mmap.const_data();
  const char* l = f + mmap.size();
  
  errno = 0;
  BatchPvalueVector tmp;
  while (errno == 0 && f && f<=(l-sizeof(tmp)) ) {
    memcpy(&tmp, f, sizeof(tmp));
    f += sizeof(tmp);
    PvalueVectorsDbRow pvecRow;
    
    pvecRow.precMz = tmp.precMz;
    pvecRow.charge = tmp.charge;
    pvecRow.scannr = tmp.scannr;
    
    pvecRow.retentionTime = tmp.retentionTime;
    pvecRow.queryCharge = tmp.queryCharge;
    
    std::vector<unsigned int> peakBins, peakScores;
    for (unsigned int j = 0; j < PvalueCalculator::kMaxScoringPeaks; ++j) {
      if (tmp.peakBins[j] != 0) {
        peakBins.push_back(tmp.peakBins[j]);
        peakScores.push_back(tmp.peakScores[j]);
      } else {
        break;
      }
    }
    
    std::vector<double> polyfit;
    for (unsigned int j = 0; j < PvalueCalculator::kPolyfitDegree + 1; ++j) {
      polyfit.push_back(tmp.polyfit[j]);
    }
    
    pvecRow.pvalCalc.initPolyfit(peakBins, peakScores, polyfit);
    
    pvalVecCollection.push_back(pvecRow);
  }
  
  if (Globals::VERB > 1) {
    std::cerr << "Read " << pvalVecCollection.size() << " pvalue vectors." << std::endl;
  }
}

void BatchPvalueVectors::parseBatchOverlapFile(
    const std::string& overlapBatchFileFN,
    std::vector< std::pair<std::string, std::string> >& overlapFNs) {
  std::ifstream pvecFileList(overlapBatchFileFN.c_str());
  std::string line;
  if (pvecFileList.is_open()) {
    while (getline(pvecFileList, line)) {
      std::string filepath1, filepath2;
      
      std::istringstream iss(line);
      iss >> filepath1 >> filepath2;
      overlapFNs.push_back(
          std::pair<std::string, std::string>(filepath1, filepath2));
    }
  }
}

void BatchPvalueVectors::processOverlapFiles(
    std::vector< std::pair<std::string, std::string> >& overlapFNs) {
  typedef std::pair<std::string, std::string> OverlapPair;
  BOOST_FOREACH(OverlapPair& p, overlapFNs) {
    std::vector<PvalueVectorsDbRow> pvalVecCollectionTail;
    std::vector<PvalueVectorsDbRow> pvalVecCollectionHead;
    
    readPvalueVectorsFile(p.first, pvalVecCollectionTail);
    readPvalueVectorsFile(p.second, pvalVecCollectionHead);
    
    batchCalculatePvaluesOverlap(pvalVecCollectionTail, pvalVecCollectionHead);
    
    remove(p.first.c_str());
    remove(p.second.c_str());
  }
  clearPvalueVectors();
}

void BatchPvalueVectors::parsePvalueVectorFile(const std::string& pvalVecInFileFN) {
  if (Globals::VERB > 2) {
    std::cerr << "Reading p-value vectors file" << std::endl;
  }
  readPvalueVectorsFile(pvalVecInFileFN, pvalVecCollection_);
  if (Globals::VERB > 2) {
    std::cerr << "Read in " << pvalVecCollection_.size() 
              << " p-value vectors from file" << std::endl;
  }
}

/* This function presumes that the pvalue vectors are sorted by precursor
   mass by writePvalueVectors() */
void BatchPvalueVectors::batchCalculatePvalues() {    
  if (Globals::VERB > 1) {
    std::cerr << "Calculating pvalues" << std::endl;
  }
  
  size_t n = pvalVecCollection_.size();
  
  time_t startTime;
  time(&startTime);
  clock_t startClock = clock();
  
#pragma omp parallel for schedule(dynamic, 1000)
  for (int i = 0; i < n; ++i) {
    if (i % 10000 == 0 && Globals::VERB > 2) {
      std::cerr << "Processing pvalue vector " << i+1 << "/" << n << " (" <<
                   i*100/n << "%)." << std::endl;
      Globals::reportProgress(startTime, startClock, i, n);
    }
    double precLimit = getUpperBound(pvalVecCollection_[i].precMz);
    std::vector<PvalueTriplet> pvalBuffer;                       
    for (size_t j = i+1; j < n; ++j) {
      if (pvalVecCollection_[j].precMz < precLimit) { 
        calculatePvalues(pvalVecCollection_[i], pvalVecCollection_[j], pvalBuffer);
      } else {
        break;
      }
    }
    pvalues_.batchWrite(pvalBuffer);
  }
  clearPvalueVectors();
  
  if (Globals::VERB > 1) {
    std::cerr << "Finished calculating pvalues." << std::endl;
  }
}

void BatchPvalueVectors::getPrecMzLimits(
    std::map<ScanId, std::pair<float, float> >& precMzLimits) {
  for (int i = 0; i < pvalVecCollection_.size(); ++i) {
    ScanId si = pvalVecCollection_[i].scannr;
    if (precMzLimits[si].first == 0.0) {
      precMzLimits[si] = std::make_pair(pvalVecCollection_[i].precMz, 
                                        pvalVecCollection_[i].precMz);
    } else if (pvalVecCollection_[i].precMz < precMzLimits[si].first) {
      precMzLimits[si].first = pvalVecCollection_[i].precMz;
    } else if (pvalVecCollection_[i].precMz > precMzLimits[si].second) {
      precMzLimits[si].second = pvalVecCollection_[i].precMz;
    }
  }
}

/* This function presumes that the pvalue vectors are sorted by precursor
   mass by writePvalueVectors() */
void BatchPvalueVectors::batchCalculateAndClusterPvalues(
    const std::string& resultTreeFN,
    const std::string& scanInfoFN) {    
  if (Globals::VERB > 1) {
    std::cerr << "Calculating pvalues" << std::endl;
  }
  
  size_t n = pvalVecCollection_.size();
  
  time_t startTime;
  time(&startTime);
  clock_t startClock = clock();
  
  std::map<ScanId, std::pair<float, float> > precMzLimits;
  if (scanInfoFN.size() > 0) {
    BatchSpectrumFiles reader;
    reader.readPrecMzLimits(scanInfoFN, precMzLimits);
  } else {
    getPrecMzLimits(precMzLimits);
  }
  
  size_t pvecBatchSize = 10000;
  size_t minPvalsForClustering = 20000000; /* = 20M */
  size_t newStartBatch = 0u, newPoisonedStartBatch = 0u;
  size_t numPvecBatches = (n - 1) / pvecBatchSize + 1;
  
  std::vector< std::vector<PvalueTriplet> > pvalBuffers(numPvecBatches);
  std::vector<bool> finishedPvalCalc(numPvecBatches);
  // MT: deque (opposed to vector) does not invalidate references!
  std::deque<ClusterJob> clusterJobs;
  ClusterJob poisonedClusterJob;
  poisonedClusterJob.finished = true;
#pragma omp parallel for schedule(dynamic)
  for (int b = 0; b < n; b += pvecBatchSize) {
    if (Globals::VERB > 2) {
      std::cerr << "Processing pvalue vector " << b+1 << "/" << n << " (" 
                << b*100/n << "%)." << std::endl;
    }
    int upperBoundIdx = (std::min)(b + pvecBatchSize, n);
    
    for (int i = b; i < upperBoundIdx; ++i) {
      double precLimit = getUpperBound(pvalVecCollection_[i].precMz);
      for (size_t j = i+1; j < n; ++j) {
        if (pvalVecCollection_[j].precMz < precLimit) { 
          calculatePvalues(pvalVecCollection_[i], pvalVecCollection_[j], 
                           pvalBuffers[b / pvecBatchSize]);
        } else {
          break;
        }
      }
    }
    finishedPvalCalc[b / pvecBatchSize] = true;
    
    bool doClustering = false;
    size_t clusterJobIdx = 0u;
  #pragma omp critical (dist_cluster)
    {
      size_t numPvals = 0u;
      size_t startIdx = newStartBatch * pvecBatchSize;
      double lowerPrecMz = pvalVecCollection_[startIdx].precMz;
      for (size_t i = newStartBatch; i < numPvecBatches; ++i) {
        if (finishedPvalCalc[i]) {
          numPvals += pvalBuffers[i].size();
          
          size_t endIdx = (std::min)((i+1) * pvecBatchSize, n) - 1;
          double upperPrecMz = pvalVecCollection_[endIdx].precMz;
          
          double threeWindows = getUpperBound(getUpperBound(getUpperBound(lowerPrecMz)));
          if (threeWindows < upperPrecMz && numPvals > minPvalsForClustering) {
            doClustering = true;
            clusterJobIdx = clusterJobs.size();
            
            ClusterJob clusterJob;
            clusterJob.startBatch = newStartBatch;
            clusterJob.endBatch = i;
            clusterJob.endIdx = endIdx;
            clusterJob.lowerPrecMz = lowerPrecMz;
            clusterJob.upperPrecMz = upperPrecMz;
            clusterJob.finished = false;
            
            newStartBatch = clusterJob.endBatch + 1;
            clusterJobs.push_back(clusterJob);
            break;
          }
        } else {
          break;
        }
      }
    }
    
    if (doClustering) {
      if (Globals::VERB > 2) {
        std::cerr << "Starting clustering job " << clusterJobIdx + 1 << std::endl;
      }
      runClusterJob(clusterJobs[clusterJobIdx], pvalBuffers, precMzLimits, resultTreeFN);
      if (Globals::VERB > 2) {
        Globals::reportProgress(startTime, startClock, clusterJobs[clusterJobIdx].endIdx, n);
      }
    }
    
    bool doPoisonedClustering = false;
  #pragma omp critical (dist_cluster)
    if (poisonedClusterJob.finished) {
      size_t numPvals = 0u;
      for (size_t i = newPoisonedStartBatch; i < clusterJobs.size(); ++i) {
        if (clusterJobs[i].finished) {
          numPvals += clusterJobs[i].poisonedPvals.size();
          if (numPvals > minPvalsForClustering && i > newPoisonedStartBatch) {
            doPoisonedClustering = true;
            
            poisonedClusterJob.startBatch = newPoisonedStartBatch;
            poisonedClusterJob.lowerPrecMz = pvalVecCollection_.front().precMz;
            poisonedClusterJob.endBatch = i;
            poisonedClusterJob.upperPrecMz = clusterJobs[i].upperPrecMz;
            poisonedClusterJob.finished = false;
            
            newPoisonedStartBatch = i;
            break;
          }
        } else {
          break;
        }
      }
    }
    
    if (doPoisonedClustering) {
      runPoisonedClusterJob(poisonedClusterJob, clusterJobs, precMzLimits, resultTreeFN);
    }
  }
  
  if (pvalVecCollection_.size() > 0) {
    std::vector<PvalueTriplet> pvalBuffer;
    if (newStartBatch < numPvecBatches) {
      for (size_t i = newStartBatch; i < numPvecBatches; ++i) {
        pvalBuffer.insert(pvalBuffer.end(), pvalBuffers[i].begin(), pvalBuffers[i].end());
        std::vector<PvalueTriplet> empty;
        pvalBuffers[i].swap(empty);
      }
      size_t startIdx = newStartBatch * pvecBatchSize;
      clusterPvals(pvalBuffer, pvalBuffer, precMzLimits, 
          pvalVecCollection_[startIdx].precMz, pvalVecCollection_.back().precMz, resultTreeFN);
    }
    
    for (size_t i = newPoisonedStartBatch; i < clusterJobs.size(); ++i) {
      pvalBuffer.insert(pvalBuffer.end(), clusterJobs[i].poisonedPvals.begin(), clusterJobs[i].poisonedPvals.end());
      std::vector<PvalueTriplet> empty;
      clusterJobs[i].poisonedPvals.swap(empty);
    }
    clusterPvals(pvalBuffer, pvalBuffer, precMzLimits, 
        pvalVecCollection_.front().precMz, pvalVecCollection_.back().precMz, resultTreeFN);
    
    pvalues_.batchWrite(pvalBuffer);
  }
  
  clearPvalueVectors();
  
  if (Globals::VERB > 1) {
    std::cerr << "Finished calculating pvalues." << std::endl;
    Globals::reportProgress(startTime, startClock, n - 1, n);
  }
}

void BatchPvalueVectors::runClusterJob(ClusterJob& clusterJob,
    std::vector< std::vector<PvalueTriplet> >& pvalBuffers,
    std::map<ScanId, std::pair<float, float> >& precMzLimits,
    const std::string& resultTreeFN) {
  std::vector<PvalueTriplet> pvalBuffer;
  for (size_t i = clusterJob.startBatch; i <= clusterJob.endBatch; ++i) {
    pvalBuffer.insert(pvalBuffer.end(), pvalBuffers[i].begin(), pvalBuffers[i].end());
    std::vector<PvalueTriplet> empty;
    pvalBuffers[i].swap(empty);
    /*if (pvalBuffer.size() > 60000000) {
      PvalueFilterAndSort::filterAndSort(pvalBuffer);
      pvalues_.batchWrite(pvalBuffer);
    }*/
  }
  
  clusterPvals(pvalBuffer, clusterJob.poisonedPvals, precMzLimits, 
      clusterJob.lowerPrecMz, clusterJob.upperPrecMz, resultTreeFN);
  clusterJob.finished = true;
  
  if (Globals::VERB > 2) {
    std::cerr << "Retained " << clusterJob.poisonedPvals.size() << " pvalues" << std::endl;
  }
}

void BatchPvalueVectors::runPoisonedClusterJob(ClusterJob& clusterJob,
    std::deque<ClusterJob>& clusterJobs,
    std::map<ScanId, std::pair<float, float> >& precMzLimits,
    const std::string& resultTreeFN) {
  std::vector<PvalueTriplet> pvalBuffer;
  for (size_t i = clusterJob.startBatch; i <= clusterJob.endBatch; ++i) {
    pvalBuffer.insert(pvalBuffer.end(), clusterJobs[i].poisonedPvals.begin(), clusterJobs[i].poisonedPvals.end());
    std::vector<PvalueTriplet> empty;
    clusterJobs[i].poisonedPvals.swap(empty);
  }
  
  std::cerr << "Starting poisoned clustering job " << clusterJob.endBatch << std::endl;
  
  clusterPvals(pvalBuffer, clusterJobs[clusterJob.endBatch].poisonedPvals, precMzLimits, 
      clusterJob.lowerPrecMz, clusterJob.upperPrecMz, resultTreeFN);
  
  if (clusterJob.startBatch == 0) {
    std::vector<PvalueTriplet> pvalBufferWrite, pvalBufferKeep;
    BOOST_FOREACH (const PvalueTriplet& pt, clusterJobs[clusterJob.endBatch].poisonedPvals) {
      if (isSafeToWrite(pt.scannr1, precMzLimits, clusterJob.upperPrecMz)
           && isSafeToWrite(pt.scannr2, precMzLimits, clusterJob.upperPrecMz)) {
        pvalBufferWrite.push_back(pt);
      } else {
        pvalBufferKeep.push_back(pt);
      }
    }
    clusterJobs[clusterJob.endBatch].poisonedPvals.swap(pvalBufferKeep);
    pvalues_.batchWrite(pvalBufferWrite);
  }
  clusterJob.finished = true;
  
  if (Globals::VERB > 2) {
    std::cerr << "Retained " << clusterJobs[clusterJob.endBatch].poisonedPvals.size() << " pvalues" << std::endl;
  }
}

void BatchPvalueVectors::clusterPvals(std::vector<PvalueTriplet>& pvalBuffer,
    std::vector<PvalueTriplet>& pvalPoisonedBuffer,
    std::map<ScanId, std::pair<float, float> >& precMzLimits, 
    float lowerPrecMz, float upperPrecMz, const std::string& resultTreeFN) {
  PvalueFilterAndSort::filter(pvalBuffer);
      
  if (Globals::VERB > 2) {
    std::cerr << "Clustering " << pvalBuffer.size() << " pvalues" << std::endl;
  }
  
  SparsePoisonedClustering matrix;
  markPoisoned(matrix, pvalBuffer, precMzLimits, lowerPrecMz, upperPrecMz);
  
  matrix.initPvals(pvalBuffer);
  matrix.setClusterPairFN(resultTreeFN);
  matrix.doClustering(dbPvalThreshold_);
  
  matrix.getPoisonedEdges(pvalPoisonedBuffer);
}

void BatchPvalueVectors::markPoisoned(SparsePoisonedClustering& matrix, 
    std::vector<PvalueTriplet>& pvalBuffer, 
    std::map<ScanId, std::pair<float, float> >& precMzLimits, 
    float lowerPrecMz, float upperPrecMz) {
  std::set<ScanId> scanIds;
  BOOST_FOREACH (const PvalueTriplet& pt, pvalBuffer) {
    scanIds.insert(pt.scannr1);
    scanIds.insert(pt.scannr2);
  }
  
  BOOST_FOREACH (const ScanId& si, scanIds) {
    if (isPoisoned(si, precMzLimits, lowerPrecMz, upperPrecMz)) {
      matrix.markPoisoned(si);
    }
  }
}

bool BatchPvalueVectors::isPoisoned(const ScanId& si,
    std::map<ScanId, std::pair<float, float> >& precMzLimits, 
    float lowerPrecMz, float upperPrecMz) {
  float minPrecMz = precMzLimits[si].first;
  float maxPrecMz = precMzLimits[si].second;
  
  return (minPrecMz < getUpperBound(lowerPrecMz)
          || upperPrecMz < getUpperBound(maxPrecMz));
}

bool BatchPvalueVectors::isSafeToWrite(const ScanId& si,
    std::map<ScanId, std::pair<float, float> >& precMzLimits, 
    float upperPrecMz) {
  float maxPrecMz = precMzLimits[si].second;
  
  return (upperPrecMz > getUpperBound(maxPrecMz));
}


/* This function presumes that the pvalue vectors are sorted by precursor
   mass by writePvalueVectors() */
void BatchPvalueVectors::batchCalculatePvaluesLibrarySearch(
    std::vector<BatchSpectrum>& querySpectra) {    
  if (Globals::VERB > 1) {
    std::cerr << "Calculating pvalues" << std::endl;
  }
  
  size_t n = pvalVecCollection_.size();
  
  time_t startTime;
  time(&startTime);
  clock_t startClock = clock();
  
  //long long numPvalsNoThresh = 0;
  int k = 0;
#pragma omp parallel for schedule(dynamic, 1000)
  for (int i = 0; i < n; ++i) {
    if (i % 10000 == 0 && Globals::VERB > 2) {
      std::cerr << "Processing pvalue vector " << i+1 << "/" << n << " (" <<
                   i*100/n << "%)." << std::endl;
      Globals::reportProgress(startTime, startClock, i, n);
    }
    std::vector<PvalueTriplet> pvalBuffer;
    double precLimitLower = getLowerBound(pvalVecCollection_[i].precMz);
    double precLimitUpper = getUpperBound(pvalVecCollection_[i].precMz);
    for (size_t j = 0; j < querySpectra.size(); ++j) {
      if (querySpectra[j].precMz < precLimitLower) {
        continue;
      }
      if (querySpectra[j].precMz < precLimitUpper) {
        calculatePvalue(pvalVecCollection_[i], querySpectra[j], pvalBuffer);
      } else {
        break;
      }
    }
    pvalues_.batchWrite(pvalBuffer);
  }
  clearPvalueVectors();
  
  if (Globals::VERB > 1) {
    //std::cerr << "Finished calculating pvalues (" << numPvalsNoThresh << " total)." << std::endl;
    std::cerr << "Finished calculating pvalues." << std::endl;
  }
  
  //PvalueFilterAndSort::filterAndSort(pvalues_.getPvaluesFN());
}

void BatchPvalueVectors::readFingerprints(
    std::vector<std::vector<unsigned short> >& mol_features, 
    std::vector<ScanId>& mol_identifiers, 
    std::vector<float>& prec_masses) {
  size_t numSpectra = pvalVecCollection_.size();
  
  unsigned int mol_count = 0;
  for (size_t i = 0; i < numSpectra; ++i) {
    PvalueVectorsDbRow s = pvalVecCollection_[i];
    
    float precMass = SpectrumHandler::calcMass(s.precMz, s.charge);
    unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(precMass);
    
    std::vector<unsigned int> peakBins = s.pvalCalc.getPeakBins();
    std::vector<unsigned short> features;
    for (unsigned int j = 0; j < (std::min)(numScoringPeaks, static_cast<unsigned int>(peakBins.size())); ++j) {
      if (peakBins[j] != 0) {
        features.push_back(peakBins[j]);
      } else {
        break;
      }
    }
    if (features.size() > 0) {
      reverse(features.begin(), features.end());
      mol_identifiers.push_back(s.scannr);
      mol_features.push_back(features);
      prec_masses.push_back(precMass);
      ++mol_count;
    }
  }
  if (Globals::VERB > 2) {
    std::cerr << "Molecules read: " << mol_count << endl;
  }
}

#ifdef FINGERPRINT_FILTER
void BatchPvalueVectors::batchCalculatePvaluesJaccardFilter() {    
  std::vector<ScanId> lib_identifiers;
  std::vector<std::vector<unsigned short> > lib_features;
  std::vector<float> lib_prec_masses;
  
  readFingerprints(lib_features, lib_identifiers, lib_prec_masses);
  
  BALL::BinaryFingerprintMethods bfm;
  float simCutoff = 0.0;
  size_t blockSize = 850;
  bfm.setCutoff(simCutoff);
  bfm.setBlockSize(blockSize);
  bfm.setVerbosityLevel(6);
  bfm.setPpmThresh(precursorTolerance_);
  bfm.setLibraryFeatures(lib_features);
  bfm.initInvertedIndicesWithPrecMasses(lib_prec_masses);
  
  if (Globals::VERB > 1) {
    std::cerr << "Calculating pvalues" << std::endl;
  }
  
  size_t n = pvalVecCollection_.size();
  size_t numBlocks = (n-1) / blockSize + 1;
  size_t numBlockPairs = numBlocks * (numBlocks+1) / 2;
  
  time_t startTime;
  time(&startTime);
  clock_t startClock = clock();
  
  //long long numPvalsNoThresh = 0;
  long long fingerPrintPairs = 0, fingerPrintComparisons = 0, pvaluePairsPassed = 0, fingerPrintPairCandidates = 0;
#pragma omp parallel for schedule(dynamic, 100)
  for (int i = 0; i < numBlockPairs; ++i) {
    int col = std::floor(std::sqrt(2 * i + 0.25) - 0.5);
    int row = i - ((col * col + col) / 2);
    if (col % 100 == 0 && row == 0 && Globals::VERB > 2) {
      std::cerr << "Processing column " << col+1 << "/" << numBlocks << " (" <<
                   col*100/numBlocks << "%)." << std::endl;
      Globals::reportProgress(startTime, startClock, col, numBlocks);
    }
    std::vector<PvalueTriplet> pvalBuffer, tmpPvalBuffer;
    bool blocksInRange = bfm.pairwiseSimilaritiesOMPThread(i, tmpPvalBuffer);
    
    int fingerPrintPairCandidatesLocal = 0;
    if (tmpPvalBuffer.size() > 0) {
      //std::cerr << i << " cand " << tmpPvalBuffer.size() << std::endl;
      
      double lowerBound = getLowerBound(pvalVecCollection_[t.scannr2.scannr].precMz);
      BOOST_FOREACH(const PvalueTriplet& t, tmpPvalBuffer) {
        if (pvalVecCollection_[t.scannr1.scannr].precMz > lowerBound) {
          calculatePvalues(pvalVecCollection_[t.scannr1.scannr], 
                           pvalVecCollection_[t.scannr2.scannr], pvalBuffer);
          ++fingerPrintPairCandidatesLocal;
        }
      }
      
      //std::cerr << i << " real " << pvalBuffer.size() << std::endl;
      pvalues_.batchWrite(pvalBuffer);
    }
    
    if (blocksInRange) {
    #pragma omp critical(count_fingerprint_pairs)
      {
        pvaluePairsPassed += pvalBuffer.size() / 2;
        fingerPrintPairs += tmpPvalBuffer.size();
        fingerPrintPairCandidates += fingerPrintPairCandidatesLocal;
        if (col != row) {
          fingerPrintComparisons += blockSize * blockSize;
        } else {
          fingerPrintComparisons += ((blockSize - 1) * blockSize) / 2;
        }
      }
    }
  }
  clearPvalueVectors();
  
  if (Globals::VERB > 1) {
    //std::cerr << "Finished calculating pvalues (" << numPvalsNoThresh << " total)." << std::endl;
    std::cerr << "Finished calculating pvalues." << std::endl;
    std::cerr << "Significant pairs found: " << pvaluePairsPassed << "/" 
              << fingerPrintPairCandidates << "/" << fingerPrintPairs 
              << "/" << fingerPrintComparisons << std::endl;
  }
  
  //PvalueFilterAndSort::filterAndSort(pvalues_.getPvaluesFN());
}
#endif

void BatchPvalueVectors::batchCalculatePvaluesOverlap(
    std::vector<PvalueVectorsDbRow>& pvalVecCollectionTail,
    std::vector<PvalueVectorsDbRow>& pvalVecCollectionHead) {
  if (Globals::VERB > 1) {
    std::cerr << "Calculating pvalues of overlap" << std::endl;
  }
  
  size_t n1 = pvalVecCollectionTail.size();
  size_t n2 = pvalVecCollectionHead.size();
  for (size_t i = 0; i < n1; ++i) {
    if (i % 10000 == 0 && Globals::VERB > 2) {
      std::cerr << "Processing pvalue vector " << i+1 << "/" << n1 << std::endl;
    }
    double precLimit = getUpperBound(pvalVecCollectionTail[i].precMz);
    std::vector<PvalueTriplet> pvalBuffer;
    for (size_t j = 0; j < n2; ++j) {
      if (pvalVecCollectionHead[j].precMz < precLimit) { 
        calculatePvalues(pvalVecCollectionTail[i], pvalVecCollectionHead[j], pvalBuffer);
      } else {
        break;
      }
    }
    
    pvalues_.batchWrite(pvalBuffer);
  }
  
  if (Globals::VERB > 1) {
    std::cerr << "Finished calculating pvalues of overlap" << std::endl;
  }
}

double BatchPvalueVectors::calculateCosineDistance(
    std::vector<unsigned int>& peakBins,
    std::vector<unsigned int>& queryPeakBins) {
  double numerator = 0.0, numerator2 = 0.0;
  std::vector<std::pair<unsigned int, double> > mziPairs, queryMziPairs;
  double div = peakBins[1];
  for (size_t i = 0; i < peakBins.size(); i += 2) {
    double normalizedIntensity = static_cast<double>(peakBins[i+1])/div;
    mziPairs.push_back(std::make_pair(peakBins[i], normalizedIntensity));
    numerator += normalizedIntensity*normalizedIntensity;
  }
  std::sort(mziPairs.begin(), mziPairs.end());
  
  div = queryPeakBins[1];
  for (size_t i = 0; i < queryPeakBins.size(); i += 2) {
    double normalizedIntensity = static_cast<double>(queryPeakBins[i+1])/div;
    queryMziPairs.push_back(std::make_pair(queryPeakBins[i], normalizedIntensity));
    numerator2 += normalizedIntensity*normalizedIntensity;
  }
  std::sort(queryMziPairs.begin(), queryMziPairs.end());
  numerator = std::sqrt(numerator) * std::sqrt(numerator2);
  
  int candIdx = 0;
  double dotProduct = 0.0;
  for (int i = 0; i < mziPairs.size(); ++i) {
    if (candIdx >= queryMziPairs.size()) break;
    while (mziPairs[i].first > queryMziPairs[candIdx].first) {
      if (++candIdx >= queryMziPairs.size()) break;
    }
    if (candIdx >= queryMziPairs.size()) break;
    if (mziPairs[i].first == queryMziPairs[candIdx].first) {
      dotProduct += mziPairs[i].second*queryMziPairs[candIdx].second;
      ++candIdx;
    }
  }
  return -100.0 * dotProduct / numerator;
}

void BatchPvalueVectors::calculatePvalues(PvalueVectorsDbRow& pvecRow, 
    PvalueVectorsDbRow& queryPvecRow, std::vector<PvalueTriplet>& pvalBuffer) {
  // skip if we are trying to score a spectrum against itself or if the charges
  // do not match
  if (queryPvecRow.scannr == pvecRow.scannr || 
      !isPvecMatch(pvecRow, queryPvecRow)) {
    return;
  }
#ifdef DOT_PRODUCT
  double cosDist = calculateCosineDistance(pvecRow.pvalCalc.getPeakBinsRef(), pvecRow.pvalCalc.getPeakBinsRef());  
  if (cosDist <= dbPvalThreshold_) {
    pvalBuffer.push_back(PvalueTriplet(std::min(pvecRow.scannr, queryPvecRow.scannr),
                                       std::max(pvecRow.scannr, queryPvecRow.scannr),
                                       cosDist));
  }
#else  
  double queryPval = queryPvecRow.pvalCalc.computePvalPolyfit(pvecRow.pvalCalc.getPeakBinsRef());
  if (queryPval <= dbPvalThreshold_) {
    double targetPval = pvecRow.pvalCalc.computePvalPolyfit(queryPvecRow.pvalCalc.getPeakBinsRef());
    if (targetPval <= dbPvalThreshold_) {
      pvalBuffer.push_back(PvalueTriplet(std::min(pvecRow.scannr, queryPvecRow.scannr),
                                         std::max(pvecRow.scannr, queryPvecRow.scannr),
                                         std::max(targetPval, queryPval)));
    }
  }
#endif
}

void BatchPvalueVectors::calculatePvalue(PvalueVectorsDbRow& pvecRow, 
                                         BatchSpectrum& querySpectrum,
                                         std::vector<PvalueTriplet>& pvalBuffer) {  
  // skip if the charges do not match
  if (pvecRow.queryCharge != querySpectrum.charge) {
    return;
  }
  
  float precMass = SpectrumHandler::calcMass(querySpectrum.precMz, 
                                             querySpectrum.charge);
  unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(precMass);
  std::vector<unsigned int> peakBins;
  peakBins.reserve(numScoringPeaks);
  for (unsigned int j = 0; j < numScoringPeaks; ++j) {
    if (querySpectrum.fragBins[j] != 0) {
      peakBins.push_back(querySpectrum.fragBins[j]);
    } else {
      break;
    }
  }
  
#ifdef DOT_PRODUCT
  double cosDist = calculateCosineDistance(pvecRow.pvalCalc.getPeakBinsRef(), 
                                           peakBins);  
  if (cosDist <= dbPvalThreshold_) {
    pvalBuffer.push_back(PvalueTriplet(pvecRow.scannr, querySpectrum.scannr,
                                       cosDist));
  }
#else
  double targetPval = pvecRow.pvalCalc.computePvalPolyfit(peakBins);
  if (targetPval <= dbPvalThreshold_) {
    pvalBuffer.push_back(PvalueTriplet(pvecRow.scannr, querySpectrum.scannr, 
                                       targetPval));
  }
#endif
}
