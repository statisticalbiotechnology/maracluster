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
#include "BatchSpectra.h"

void BatchPvalueVectors::initPvecRow(const MassChargeCandidate& mcc, 
                                    const BatchSpectrum& spec,
                                    PvalueVectorsDbRow& pvecRow) {
  pvecRow.precMass = mcc.mass;
  pvecRow.charge = mcc.charge;
  pvecRow.scannr = spec.scannr;
  pvecRow.retentionTime = spec.retentionTime;
  
  unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(spec.precMass);
#ifdef DOT_PRODUCT
  numScoringPeaks *= 2;
  pvecRow.peakBins.reserve(numScoringPeaks);
  for (unsigned int j = 0; j < numScoringPeaks; ++j) {
    if (spec.fragBins[j] != 0) {
      pvecRow.peakBins.push_back(spec.fragBins[j]);
    } else if (j % 2 == 0) {
      break;
    }
  }
#else
  pvecRow.peakBins.reserve(numScoringPeaks);
  for (unsigned int j = 0; j < numScoringPeaks; ++j) {
    if (spec.fragBins[j] != 0) {
      pvecRow.peakBins.push_back(spec.fragBins[j]);
    } else {
      break;
    }
  }
#endif
}

void BatchPvalueVectors::insertMassChargeCandidate(
    MassChargeCandidate& mcc, BatchSpectrum& spec) {
  if (BatchGlobals::VERB > 4) {
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
  
  if (BatchGlobals::VERB > 4) {
    std::cerr << "Inserted mass charge candidate into database" << std::endl;
  }
}

void BatchPvalueVectors::batchInsert(PeakCounts& peakCounts, bool forceInsert) {
  if (pvalVecBatch_.size() >= 5000 || forceInsert) {
    if (BatchGlobals::VERB > 3) {
      std::cerr << "Inserting " << pvalVecBatch_.size() << " spectra into database" << std::endl;
    }
    //std::cerr << "Calculating " << pvalVecBatch_.size() << " p-value vectors" << std::endl;
  #pragma omp parallel for schedule(dynamic, 100)
    for (size_t i = 0; i < pvalVecBatch_.size(); ++i) {
      calculatePvalueVector(pvalVecBatch_[i], peakCounts);
    }
    pvalVecBatch_.clear();
    
    if (BatchGlobals::VERB > 3) {
      std::cerr << "Currently there are " << pvalVecCollection_.size() << " spectra in the database" << std::endl;
    }
  }
}

void BatchPvalueVectors::initPvalCalc(PvalueCalculator& pvalCalc, 
    PvalueVectorsDbRow& pvecRow, PeakCounts& peakCounts, 
    const int numQueryPeaks, const bool polyfit) {
#ifdef DOT_PRODUCT
  pvalCalc.init(pvecRow.peakBins, std::vector<double>() );
#else
  PeakDistribution distribution;
  double precMz = SpectrumHandler::calcPrecMz(pvecRow.precMass, pvecRow.queryCharge);
  peakCounts.generatePeakDistribution(precMz, pvecRow.queryCharge, 
                                      distribution, numQueryPeaks);
  
  pvalCalc.initFromPeakBins(pvecRow.peakBins, distribution.getDistribution());
  
  if (polyfit) {
    pvalCalc.computePvalVectorPolyfit();
  } else {
    pvalCalc.computePvalVector();
  } 
#endif
}

void BatchPvalueVectors::calculatePvalueVector(PvalueVectorsDbRow& pvecRow,
    PeakCounts& peakCounts) {
  if (BatchGlobals::VERB > 4) {
    std::cerr << "Inserting pvalue vector " << pvecRow.scannr << std::endl;
  }
  bool polyfit = true;
  unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(pvecRow.precMass);
  initPvalCalc(pvecRow.pvalCalc, pvecRow, peakCounts, numScoringPeaks, polyfit);
  
  if (pvecRow.pvalCalc.getNumScoringPeaks() >= PvalueCalculator::getMinScoringPeaks(pvecRow.precMass)) {    
  #pragma omp critical (store_pvec)
    {
      pvalVecCollection_.push_back(pvecRow);
    }
        
    if (BatchGlobals::VERB > 4) {
      std::cerr << "Inserted pvalue vector " << pvecRow.scannr << std::endl;
    }
  } else if (BatchGlobals::VERB > 2) {
    std::cerr << "Ignoring " << pvecRow.scannr << ". Not enough scoring peaks: " << pvecRow.pvalCalc.getNumScoringPeaks() << std::endl;
  }
}

void BatchPvalueVectors::sortPvalueVectors() {
  if (BatchGlobals::VERB > 2) {
    std::cerr << "Sorting pvalue vectors" << std::endl;
  }
  std::sort(pvalVecCollection_.begin(), pvalVecCollection_.end());
}

void BatchPvalueVectors::writePvalueVectors(
    const std::string& pvalueVectorsBaseFN) {
  std::string pvalueVectorsFN = pvalueVectorsBaseFN + ".dat";
  std::string pvalueVectorsHeadFN = pvalueVectorsBaseFN + ".head.dat";
  std::string pvalueVectorsTailFN = pvalueVectorsBaseFN + ".tail.dat";
  
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Writing pvalue vectors" << std::endl;
  }

  std::vector<BatchPvalueVector> headList, tailList, allList;
  double headOverlapLimit = getLowerBound(pvalVecCollection_.front().precMass);
  double tailOverlapLimit = getUpperBound(pvalVecCollection_.back().precMass);
  size_t n = pvalVecCollection_.size();
  for (size_t i = 0; i < n; ++i) {
    if (i % 100000 == 0 && BatchGlobals::VERB > 2) {
      std::cerr << "Writing pvalue vector " << i << "/" << n << std::endl;
    }

    insert(pvalVecCollection_[i], allList);
    if (pvalVecCollection_[i].precMass < headOverlapLimit) {
      insert(pvalVecCollection_[i], headList);
    }
    if (pvalVecCollection_[i].precMass > tailOverlapLimit) {
      insert(pvalVecCollection_[i], tailList);
    }
  }

  bool append = false;
  BinaryInterface::write<BatchPvalueVector>(allList, pvalueVectorsFN, append);
  BinaryInterface::write<BatchPvalueVector>(headList, pvalueVectorsHeadFN, append);
  BinaryInterface::write<BatchPvalueVector>(tailList, pvalueVectorsTailFN, append);
  
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Finished writing pvalue vectors" << std::endl;
  }
}

void BatchPvalueVectors::insert(PvalueVectorsDbRow& pvecRow, std::vector<BatchPvalueVector>& pvecList) {
  if (BatchGlobals::VERB > 4) {
    std::cerr << "Inserting pvalue vector into pvalue vectors " <<
                 "table asynchronously" << std::endl;
  }  
  BatchPvalueVector pvec;
  
  pvec.precMass = pvecRow.precMass;
  pvec.charge = pvecRow.charge;
  pvec.scannr = pvecRow.scannr;
  
  pvec.retentionTime = pvecRow.retentionTime;
  pvec.queryCharge = pvecRow.queryCharge;
  
  pvecRow.pvalCalc.copyPolyfit(pvec.peakBins, pvec.peakScores, pvec.polyfit);
  
  pvecList.push_back(pvec);
  
  if (BatchGlobals::VERB > 4) {
    std::cerr << "Put pvalue vector insertion into queue" << std::endl;
  }
}

void BatchPvalueVectors::readPvalueVectorsFile(const std::string& pvalueVectorsFN,
    std::vector<PvalueVectorsDbRow>& pvalVecCollection) {  
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Reading in pvalue vectors from " << pvalueVectorsFN << std::endl;
  }

  if (!BatchGlobals::fileExists(pvalueVectorsFN)) {
    std::cerr << "Ignoring missing file " << pvalueVectorsFN << std::endl;
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
    
    pvecRow.precMass = tmp.precMass;
    pvecRow.charge = tmp.charge;
    pvecRow.scannr = tmp.scannr;
    
    pvecRow.retentionTime = tmp.retentionTime;
    pvecRow.queryCharge = tmp.queryCharge;
    
    std::vector<unsigned int> peakScores;
    for (unsigned int j = 0; j < PvalueCalculator::kMaxScoringPeaks; ++j) {
      if (tmp.peakBins[j] != 0) {
        pvecRow.peakBins.push_back(tmp.peakBins[j]);
        peakScores.push_back(tmp.peakScores[j]);
      } else {
        break;
      }
    }
    
    std::vector<double> polyfit;
    for (unsigned int j = 0; j < PvalueCalculator::kPolyfitDegree + 1; ++j) {
      polyfit.push_back(tmp.polyfit[j]);
    }
    
    pvecRow.pvalCalc.initPolyfit(pvecRow.peakBins, peakScores, polyfit);
    
    pvalVecCollection.push_back(pvecRow);
  }
  
  if (BatchGlobals::VERB > 1) {
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
  }
  clearPvalueVectors();
  
  //PvalueFilterAndSort::filterAndSort(pvalues_.getPvaluesFN());
}

void BatchPvalueVectors::parsePvalueVectorFile(const std::string& pvalVecInFileFN) {
  if (BatchGlobals::VERB > 2) {
    std::cerr << "Reading p-value vectors file" << std::endl;
  }
  readPvalueVectorsFile(pvalVecInFileFN, pvalVecCollection_);
  if (BatchGlobals::VERB > 2) {
    std::cerr << "Read in " << pvalVecCollection_.size() << " p-value vectors from file" << std::endl;
  }
}

/* This function presumes that the pvalue vectors are sorted by precursor
   mass by writePvalueVectors() */
void BatchPvalueVectors::batchCalculatePvalues() {    
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Calculating pvalues" << std::endl;
  }
  
  size_t n = pvalVecCollection_.size();
  
  time_t startTime;
  time(&startTime);
  clock_t startClock = clock();
  
  //long long numPvalsNoThresh = 0;
#pragma omp parallel for schedule(dynamic, 1000)
  for (size_t i = 0; i < n; ++i) {
    if (i % 10000 == 0 && BatchGlobals::VERB > 2) {
      std::cerr << "Processing pvalue vector " << i+1 << "/" << n << " (" <<
                   i*100/n << "%)." << std::endl;
      BatchGlobals::reportProgress(startTime, startClock, i, n);
    }
    double precLimit = getUpperBound(pvalVecCollection_[i].precMass);
    std::vector<PvalueTriplet> pvalBuffer;                       
    for (size_t j = i+1; j < n; ++j) {
      if (pvalVecCollection_[j].precMass < precLimit) { 
        calculatePvalues(pvalVecCollection_[i], pvalVecCollection_[j], pvalBuffer);
        //#pragma omp atomic
        //++numPvalsNoThresh;
      } else {
        break;
      }
    }
    pvalues_.batchWrite(pvalBuffer);
  }
  clearPvalueVectors();
  
  if (BatchGlobals::VERB > 1) {
    //std::cerr << "Finished calculating pvalues (" << numPvalsNoThresh << " total)." << std::endl;
    std::cerr << "Finished calculating pvalues." << std::endl;
  }
  
  //PvalueFilterAndSort::filterAndSort(pvalues_.getPvaluesFN());
}

/* This function presumes that the pvalue vectors are sorted by precursor
   mass by writePvalueVectors() */
void BatchPvalueVectors::batchCalculatePvaluesLibrarySearch(
    std::vector<BatchSpectrum>& querySpectra) {    
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Calculating pvalues" << std::endl;
  }
  
  size_t n = pvalVecCollection_.size();
  
  time_t startTime;
  time(&startTime);
  clock_t startClock = clock();
  
  //long long numPvalsNoThresh = 0;
  int k = 0;
#pragma omp parallel for schedule(dynamic, 1000)
  for (size_t i = 0; i < n; ++i) {
    if (i % 10000 == 0 && BatchGlobals::VERB > 2) {
      std::cerr << "Processing pvalue vector " << i+1 << "/" << n << " (" <<
                   i*100/n << "%)." << std::endl;
      BatchGlobals::reportProgress(startTime, startClock, i, n);
    }
    std::vector<PvalueTriplet> pvalBuffer;
    double precLimitLower = getLowerBound(pvalVecCollection_[i].precMass);
    double precLimitUpper = getUpperBound(pvalVecCollection_[i].precMass);
    for (size_t j = 0; j < querySpectra.size(); ++j) {
      if (querySpectra[j].precMass < precLimitLower) {
        continue;
      }
      if (querySpectra[j].precMass < precLimitUpper) { 
        calculatePvalue(pvalVecCollection_[i], querySpectra[j], pvalBuffer);
      } else {
        break;
      }
    }
    pvalues_.batchWrite(pvalBuffer);
  }
  clearPvalueVectors();
  
  if (BatchGlobals::VERB > 1) {
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
    
    unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(s.precMass);
    
    std::vector<unsigned short> features;
    for (unsigned int j = 0; j < (std::min)(numScoringPeaks, static_cast<unsigned int>(s.peakBins.size())); ++j) {
      if (s.peakBins[j] != 0) {
        features.push_back(s.peakBins[j]);
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
  if (BatchGlobals::VERB > 2) {
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
  
  if (BatchGlobals::VERB > 1) {
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
  for (size_t i = 0; i < numBlockPairs; ++i) {
    int col = std::floor(std::sqrt(2 * i + 0.25) - 0.5);
    int row = i - ((col * col + col) / 2);
    if (col % 100 == 0 && row == 0 && BatchGlobals::VERB > 2) {
      std::cerr << "Processing column " << col+1 << "/" << numBlocks << " (" <<
                   col*100/numBlocks << "%)." << std::endl;
      BatchGlobals::reportProgress(startTime, startClock, col, numBlocks);
    }
    std::vector<PvalueTriplet> pvalBuffer, tmpPvalBuffer;
    bool blocksInRange = bfm.pairwiseSimilaritiesOMPThread(i, tmpPvalBuffer);
    
    int fingerPrintPairCandidatesLocal = 0;
    if (tmpPvalBuffer.size() > 0) {
      //std::cerr << i << " cand " << tmpPvalBuffer.size() << std::endl;
      
      double lowerBound = getLowerBound(pvalVecCollection_[t.scannr2.scannr].precMass);
      BOOST_FOREACH(const PvalueTriplet& t, tmpPvalBuffer) {
        if (pvalVecCollection_[t.scannr1.scannr].precMass > lowerBound) {
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
  
  if (BatchGlobals::VERB > 1) {
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
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Calculating pvalues of overlap" << std::endl;
  }
  
  size_t n1 = pvalVecCollectionTail.size();
  size_t n2 = pvalVecCollectionHead.size();
  for (size_t i = 0; i < n1; ++i) {
    if (i % 10000 == 0 && BatchGlobals::VERB > 2) {
      std::cerr << "Processing pvalue vector " << i+1 << "/" << n1 << std::endl;
    }
    double precLimit = getUpperBound(pvalVecCollectionTail[i].precMass);
    std::vector<PvalueTriplet> pvalBuffer;
    for (size_t j = 0; j < n2; ++j) {
      //std::cerr << pvalVecCollectionHead[j].precMass << " " << precLimit << std::endl;
      if (pvalVecCollectionHead[j].precMass < precLimit) { 
        //std::cerr << pvalVecCollectionHead[j].precMass << " " << precLimit << std::endl;
        calculatePvalues(pvalVecCollectionTail[i], pvalVecCollectionHead[j], pvalBuffer);
      } else {
        break;
      }
    }
    
    pvalues_.batchWrite(pvalBuffer);
  }
  
  if (BatchGlobals::VERB > 1) {
    std::cerr << "Finished calculating pvalues of overlap" << std::endl;
  }
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
  double numerator = 0.0, numerator2 = 0.0;
  std::vector<std::pair<unsigned int, double> > mziPairs, queryMziPairs;
  double div = pvecRow.peakBins[1];
  for (size_t i = 0; i < pvecRow.peakBins.size(); i += 2) {
    double normalizedIntensity = static_cast<double>(pvecRow.peakBins[i+1])/div;
    mziPairs.push_back(std::make_pair(pvecRow.peakBins[i], normalizedIntensity));
    numerator += normalizedIntensity*normalizedIntensity;
  }
  std::sort(mziPairs.begin(), mziPairs.end());
  
  div = queryPvecRow.peakBins[1];
  for (size_t i = 0; i < queryPvecRow.peakBins.size(); i += 2) {
    double normalizedIntensity = static_cast<double>(queryPvecRow.peakBins[i+1])/div;
    queryMziPairs.push_back(std::make_pair(queryPvecRow.peakBins[i], normalizedIntensity));
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
  double cosDist = -100.0 * dotProduct / numerator;
  if (cosDist < dbPvalThreshold_) {
    pvalBuffer.push_back(PvalueTriplet(std::min(pvecRow.scannr, queryPvecRow.scannr),
                                       std::max(pvecRow.scannr, queryPvecRow.scannr),
                                       cosDist));
    //pvalBuffer.push_back(PvalueTriplet(queryPvecRow.scannr, pvecRow.scannr, cosDist));
  }
#else  
  double queryPval = queryPvecRow.pvalCalc.computePvalPolyfit(pvecRow.peakBins);
  if (queryPval < dbPvalThreshold_) {
    double targetPval = pvecRow.pvalCalc.computePvalPolyfit(queryPvecRow.peakBins);
    if (targetPval < dbPvalThreshold_) {
      pvalBuffer.push_back(PvalueTriplet(std::min(pvecRow.scannr, queryPvecRow.scannr),
                                         std::max(pvecRow.scannr, queryPvecRow.scannr),
                                         std::max(targetPval, queryPval)));
      //pvalBuffer.push_back(PvalueTriplet(queryPvecRow.scannr, pvecRow.scannr, queryPval));
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
  
  unsigned int numScoringPeaks = PvalueCalculator::getMaxScoringPeaks(querySpectrum.precMass);
  std::vector<unsigned int> peakBins;
  peakBins.reserve(numScoringPeaks);
  for (unsigned int j = 0; j < numScoringPeaks; ++j) {
    if (querySpectrum.fragBins[j] != 0) {
      peakBins.push_back(querySpectrum.fragBins[j]);
    } else {
      break;
    }
  }
  
  double targetPval = pvecRow.pvalCalc.computePvalPolyfit(peakBins);
  if (targetPval < dbPvalThreshold_) {
    pvalBuffer.push_back(PvalueTriplet(pvecRow.scannr, querySpectrum.scannr, 
                                       targetPval));
  }
}
