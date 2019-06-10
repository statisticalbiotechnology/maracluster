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

#include "MSFileMerger.h"
#include "Version.h"

namespace maracluster {

using pwiz::msdata::MSDataPtr;
using pwiz::msdata::MSData;
using pwiz::msdata::MSDataFile;
using pwiz::msdata::MSDataMerger;
using pwiz::msdata::SpectrumListSimplePtr;
using pwiz::msdata::SpectrumListSimple;
using pwiz::msdata::SpectrumListPtr;
using pwiz::msdata::SpectrumPtr;

using pwiz::msdata::Software;
using pwiz::msdata::SoftwarePtr;
using pwiz::msdata::SourceFile;
using pwiz::msdata::SourceFilePtr;
using pwiz::msdata::DataProcessing;
using pwiz::msdata::DataProcessingPtr;
using pwiz::msdata::ProcessingMethod;

int MSFileMerger::maxMSFilePtrs_ = 40; // memory constrained
int MSFileMerger::maxSpectraPerFile_ = 200000; // memory constrained
unsigned int MSFileMerger::maxConsensusSpectraPerFile_ = 200000u; // search engine memory constrained

int MSFileMerger::mergeMethod_ = 10;
bool MSFileMerger::normalize_ = true;

void MSFileMerger::parseClusterFileForMerge(const std::string& clusterFile,
    const size_t minClusterSize) {
  std::ifstream clusterStream(clusterFile.c_str());
  std::string line;
  if (clusterStream.is_open()) {
    bool newCluster = true;
    size_t mergeCount = 0;
    while (getline(clusterStream, line)) {
      std::istringstream iss(line);
      std::string filepath;
      unsigned int scannr = 0;

      iss >> filepath >> scannr;

      if (iss.fail()) {
        newCluster = true;
        iss.clear();
      } else {
        if (newCluster) {
          ScanMergeInfoSet s;
          s.mergedScanId = fileList_.getScanId(spectrumOutFN_, ++mergeCount);
          if (combineSets_.size() > 0 && combineSets_.back().size() < minClusterSize) {
            combineSets_.pop_back();
          }
          combineSets_.push_back(s);
          newCluster = false;
        }

        ScanId globalIdx = fileList_.getScanId(filepath, scannr);
        combineSets_.back().push_back(globalIdx);
      }
    }

    if (combineSets_.size() > 0 && combineSets_.back().size() < minClusterSize) {
      combineSets_.pop_back();
    }
  }
}

void MSFileMerger::mergeAllSpectra(const std::string& spectrumInFN) {
  MSReaderList readerList;
  MSDataFile msd(spectrumInFN, &readerList);
  SpectrumListPtr sl = msd.run.spectrumListPtr;

  MSClusterMerge::init();

  std::vector<SpectrumPtr> spectra;
  for (unsigned int j = 0; j < sl->size(); ++j) {
    SpectrumPtr s = sl->spectrum(j, true);
    spectra.push_back(s);
  }

  SpectrumListSimplePtr mergedSpectra(new SpectrumListSimple);
  mergedSpectra->dp = DataProcessingPtr(new DataProcessing("MaRaCluster_consensus_builder"));

  ScanId mergedScanId = fileList_.getScanId(spectrumOutFN_,1u);
  mergeSpectraSetMSCluster(spectra, mergedScanId, mergedSpectra);

  std::cerr << "Creating merged MSData file" << std::endl;
  MSData msdMerged;
  msdMerged.id = msdMerged.run.id = "merged_spectra";
  msdMerged.run.spectrumListPtr = mergedSpectra;

  writeMSData(msdMerged, spectrumOutFN_);
}

void MSFileMerger::mergeSpectraSetMSCluster(
    std::vector<SpectrumPtr>& spectra,
    ScanId scanId, SpectrumListSimplePtr mergedSpectra) {
  std::vector< std::vector<BinnedMZIntensityPair> > cluster;
  std::vector<MassChargeCandidate> allMccs;
  bool first = true;
  SpectrumPtr consensusSpec;
  BOOST_FOREACH(SpectrumPtr s, spectra) {
    if (first) {
      consensusSpec = SpectrumPtr(new pwiz::msdata::Spectrum(*s));
      first = false;
    }

    std::vector<MassChargeCandidate> mccs;
    SpectrumHandler::getMassChargeCandidates(s, mccs);
    allMccs.insert(allMccs.end(), mccs.begin(), mccs.end());

    std::vector<MZIntensityPair> mziPairs;
    SpectrumHandler::getMZIntensityPairs(s, mziPairs);
    if (normalize_) SpectrumHandler::normalizeIntensitiesMSCluster(mziPairs);

    std::vector<BinnedMZIntensityPair> binnedMziPairs;
    MSClusterMerge::binMZIntensityPairs(mziPairs, binnedMziPairs);
    cluster.push_back(binnedMziPairs);
  }
  std::vector<BinnedMZIntensityPair> mergedMziPairs;
  MSClusterMerge::merge(cluster, mergedMziPairs);

  std::vector<MZIntensityPair> mziPairs;
  MSClusterMerge::unbinMZIntensityPairs(mergedMziPairs, mziPairs);

  SpectrumHandler::setMZIntensityPairs(consensusSpec, mziPairs);

  SpectrumHandler::fixMetaData(consensusSpec);
  consensusSpec->scanList.set(pwiz::cv::MS_mean_of_spectra);

  std::vector<MassChargeCandidate> consensusMccs;
  mergeMccs(allMccs, consensusMccs, scanId.scannr);
  /*consensusMccs = allMccs;*/

  addSpectrumWithMccs(consensusSpec, consensusMccs, scanId.scannr, mergedSpectra);
}

void MSFileMerger::mergeMccs(std::vector<MassChargeCandidate>& allMccs,
      std::vector<MassChargeCandidate>& consensusMccs, int clusterIdx) {
  MSClusterMerge::mergeMccs(allMccs, consensusMccs);
}

void MSFileMerger::mergeSpectra() {
  if (fileList_.getFilePaths().size() > 0 /*maxMSFilePtrs_*/) {
    mergeSpectraScalable();
  } else {
    mergeSpectraSmall();
  }
}

/* Merges spectra in-memory */
void MSFileMerger::mergeSpectraSmall() {
  /* create the MSData object in memory */
  std::cerr << "Creating container for merged spectra" << std::endl;
  MSData msdMerged;
  msdMerged.id = msdMerged.run.id = "merged_spectra";

  SpectrumListSimplePtr mergedSpectra(new SpectrumListSimple);
  mergedSpectra->dp = DataProcessingPtr(new DataProcessing("MaRaCluster_consensus_builder"));

  std::cerr << "Merging spectra" << std::endl;
  int count = 0;
  int partIdx = 0;

  std::map<std::string, MSDataPtr> msDataPtrs;
  MSClusterMerge::init();

  std::sort(combineSets_.begin(), combineSets_.end(), ScanMergeInfoSet::lowerScannr);

  BOOST_FOREACH (ScanMergeInfoSet& mergeSet, combineSets_) {
    std::vector<SpectrumPtr> spectra;
    BOOST_FOREACH (ScanId scanId, mergeSet.scans) {
      std::string filepath = fileList_.getFilePath(scanId);
      unsigned int scannr = fileList_.getScannr(scanId);

      if (msDataPtrs.find(filepath) == msDataPtrs.end()) {
        MSDataPtr msd(new MSDataFile(filepath));
        msDataPtrs[filepath] = msd;
        std::cerr << "Extracting spectra from " << filepath
                  << ", numFilePointers:" << msDataPtrs.size() << std::endl;
      }

      SpectrumListPtr sl = msDataPtrs[filepath]->run.spectrumListPtr;

      size_t result = getSpectrumIdxFromScannr(sl, scannr);
      SpectrumPtr s = sl->spectrum(result, true);
      spectra.push_back(s);
    }
    mergeSpectraSetMSCluster(spectra, mergeSet.mergedScanId, mergedSpectra);

    if (++count % 5000 == 0)
      std::cerr << "Merged " << count << "/" << combineSets_.size() << std::endl;

    if (mergedSpectra->spectra.size() >= maxConsensusSpectraPerFile_
          || count == static_cast<int>(combineSets_.size())) {
      std::cerr << "Creating merged MSData file" << std::endl;
      msdMerged.run.spectrumListPtr = mergedSpectra;

      std::cerr << "#spectra: " << msdMerged.run.spectrumListPtr->size() << std::endl;
      ++partIdx;

      typedef std::pair<const std::string, MSDataPtr> FnMsDatPtrPair;
      BOOST_FOREACH (FnMsDatPtrPair& fnMsDatPtrPair, msDataPtrs) {
        msdMerged.fileDescription.sourceFilePtrs.push_back(
            fnMsDatPtrPair.second->fileDescription.sourceFilePtrs.front());
      }

      std::string partSpecOutFN = getPartFN(
          spectrumOutFN_, "part" + boost::lexical_cast<std::string>(partIdx));
      writeMSData(msdMerged, partSpecOutFN);

      mergedSpectra = SpectrumListSimplePtr(new SpectrumListSimple);
      mergedSpectra->dp = DataProcessingPtr(new DataProcessing("MaRaCluster_consensus_builder"));
    }
  }
}


/**
 * Terminology:
 *   Batch: determines how many input files are processed before writing intermediate files
 *          (depends on maxMSFilePtrs_)
 *   Bin: determines in how many files the intermediate files are split
 *        (depends on maxSpectraPerFile_ and maxConsensusSpectraPerFile_)
 **/
void MSFileMerger::mergeSpectraScalable() {
  std::map<ScanId, ScanId> scannrToMergedScannr;
  createScannrToMergedScannrMap(scannrToMergedScannr);

  numBatches_ = calcNumBatches(fileList_.size(), maxMSFilePtrs_);
  numMSFilePtrsPerBatch_ = calcNumBatches(fileList_.size(), numBatches_);
  
  size_t numClusterBinsIn = calcNumBatches(
      scannrToMergedScannr.size(), maxSpectraPerFile_);
  size_t numClusterBinsOut = calcNumBatches(
      combineSets_.size(), maxConsensusSpectraPerFile_);
  numClusterBins_ = std::max(numClusterBinsIn, numClusterBinsOut);

  splitSpecFilesByConsensusSpec(scannrToMergedScannr);
  mergeSplitSpecFiles();
}

void MSFileMerger::createScannrToMergedScannrMap(
    std::map<ScanId, ScanId>& scannrToMergedScannr) {
  BOOST_FOREACH (ScanMergeInfoSet& mergeSet, combineSets_) {
    BOOST_FOREACH (ScanId scannr, mergeSet.scans) {
      scannrToMergedScannr[scannr] = mergeSet.mergedScanId;
    }
  }
}

void MSFileMerger::splitSpecFilesByConsensusSpec(
    std::map<ScanId, ScanId>& scannrToMergedScannr) {
  for (size_t batchNr = 0; batchNr < numBatches_; ++batchNr) {
    std::vector<SpectrumListSimplePtr> spectrumListsAcc(numClusterBins_);
    for (size_t k = 0; k < numClusterBins_; ++k) {
      spectrumListsAcc[k] = SpectrumListSimplePtr(new SpectrumListSimple);
      spectrumListsAcc[k]->dp = DataProcessingPtr(new DataProcessing("MaRaCluster_consensus_builder"));
    }
    int startFileIdx = batchNr*numMSFilePtrsPerBatch_;
    int endFileIdx = (std::min)((batchNr+1)*numMSFilePtrsPerBatch_, fileList_.size());
    std::vector<MSDataPtr> msDataPtrs(endFileIdx - startFileIdx);
    int removeIdx = -1; /* index for the output mergeFile */
  #pragma omp parallel 
    {
      /* create one spectrumLists 2D vector for each thread */
      std::vector<std::vector<SpectrumPtr> > spectrumLists(numClusterBins_);
      
    #pragma omp for schedule(dynamic, 1)
      for (int i = startFileIdx; i < endFileIdx; ++i) {
        std::string filePath = fileList_.getFilePath(i);
        if (filePath == spectrumOutFN_) {
          removeIdx = i - startFileIdx;
          continue;
        }
        std::cerr << "Splitting " << filePath
                  << " (" << i*100 / (fileList_.size()-1) << "%)" << std::endl;
        
        MSReaderList readerList;
      #pragma omp critical (load_msdata_file)
        {
          /* loading multiple MSDataFile objects in parallel is slower (15-20%) 
             and (in rare cases) seems to cause runtime exceptions */
          msDataPtrs.at(i - startFileIdx) = MSDataPtr(new MSDataFile(filePath, &readerList));
        }
        
        SpectrumListPtr sl = msDataPtrs.at(i - startFileIdx)->run.spectrumListPtr;
          
        /* use a skeleton (i.e. without spectra data) in the msDataPtrs vector 
           to propagate meta data to the final file */
        msDataPtrs.at(i - startFileIdx)->run.spectrumListPtr = SpectrumListSimplePtr(new SpectrumListSimple);
        
        for (unsigned int j = 0; j < sl->size(); ++j) {
          SpectrumPtr s = sl->spectrum(j, true);
          unsigned int scannr = SpectrumHandler::getScannr(s);
          ScanId scanId = fileList_.getScanId(filePath, scannr);
          if (scannrToMergedScannr.find(scanId) != scannrToMergedScannr.end()) {
            SpectrumHandler::setScannr(s, scanId);
            SpectrumHandler::updateMassChargeCandidates(s); // fixes some incompatibility issues between mgf and ms2

            ScanId mergedScannr = scannrToMergedScannr[scanId];
            unsigned int clusterBin = getClusterBin(mergedScannr);
            spectrumLists[clusterBin].push_back(s);
          }
        }
      }
      
      /* collect spectrumLists from all threads */
    #pragma omp critical (merge_speclists)
      {
        for (size_t k = 0; k < numClusterBins_; ++k) {
          spectrumListsAcc[k]->spectra.reserve(
              spectrumListsAcc[k]->spectra.size() + spectrumLists[k].size());
          spectrumListsAcc[k]->spectra.insert(
              spectrumListsAcc[k]->spectra.end(),
              spectrumLists[k].begin(),
              spectrumLists[k].end() );
        }
      }
    } /* end omp parallel */
     
    if (removeIdx >= 0) {
      msDataPtrs.erase(msDataPtrs.begin() + removeIdx);
    }
    
    writeClusterBins(batchNr, msDataPtrs, spectrumListsAcc);
  }
  std::cerr << "Finished splitting ms2 files" << std::endl;
}

void MSFileMerger::mergeSplitSpecFiles() {
  MSClusterMerge::init();

  SpectrumListSimplePtr mergedSpectra(new SpectrumListSimple);

  SoftwarePtr softwarePtr = SoftwarePtr(new Software("MaRaCluster"));
  softwarePtr->version = VERSION;
  // FIXME: get proper PSI:MS Id for MaRaCluster?
  softwarePtr->set(pwiz::cv::MS_custom_unreleased_software_tool);
  DataProcessingPtr dpPtr = DataProcessingPtr(new DataProcessing("MaRaCluster_consensus_builder"));
  ProcessingMethod pm;
  pm.softwarePtr = softwarePtr;
  pm.set(pwiz::cv::MS_precursor_recalculation);
  dpPtr->processingMethods.push_back(pm);
  mergedSpectra->dp = dpPtr;

  unsigned int partIdx = 0u;
  
  std::vector<MSDataPtr> msdVector;
  
  for (size_t clusterBin = 0; clusterBin < numClusterBins_; ++clusterBin) {
    std::cerr << "Merging spectra in bin " << clusterBin+1 << "/" << numClusterBins_ << std::endl;
    
    mergeSpectraBin(clusterBin, msdVector, mergedSpectra);

    while (mergedSpectra->size() > maxConsensusSpectraPerFile_
        || (mergedSpectra->size() > 0u && clusterBin == numClusterBins_ - 1)) {
      std::cerr << "Creating merged MSData file" << std::endl;
      MSDataMerger msdMerged(msdVector);
      msdMerged.id = msdMerged.run.id = "consensus_spectra";
      msdMerged.softwarePtrs.push_back(softwarePtr);

      SpectrumListSimplePtr writeSpectra(new SpectrumListSimple);
      writeSpectra->dp = dpPtr;
      
      std::sort(mergedSpectra->spectra.begin(), mergedSpectra->spectra.end(), SpectrumHandler::lessScannr);
      
      size_t idx = 0;
      BOOST_FOREACH (SpectrumPtr& s, mergedSpectra->spectra) {
        s->index = idx;
        writeSpectra->spectra.push_back(s);
        if (++idx >= maxConsensusSpectraPerFile_) break;
      }

      msdMerged.run.spectrumListPtr = writeSpectra;

      std::string partSpecOutFN = getPartFN(spectrumOutFN_, "part" +
          boost::lexical_cast<std::string>(++partIdx));
      writeMSData(msdMerged, partSpecOutFN);

      mergedSpectra->spectra.erase(mergedSpectra->spectra.begin(), mergedSpectra->spectra.begin() + idx);
    }
    msdVector.clear();
  }
}

void MSFileMerger::mergeSpectraBin(size_t clusterBin, std::vector<MSDataPtr>& msdVector,
    SpectrumListSimplePtr mergedSpectra) {
  {
    for (unsigned int i = 0; i < numBatches_; ++i) {
      std::string partSpecOutFN = getPartFN(spectrumOutFN_, "part" +
            boost::lexical_cast<std::string>(clusterBin) + "_" +
            boost::lexical_cast<std::string>(i) + "_tmpfile");
      MSDataPtr msd(new MSDataFile(partSpecOutFN));
      msdVector.push_back(msd);
    }

    std::vector< std::vector<MergeScanIndex> > scanIndices(numBatches_);
    std::vector< std::pair< std::vector<SpectrumPtr>, ScanId> > spectra;

    for (unsigned int i = 0; i < combineSets_.size(); ++i) {
      if (getClusterBin(combineSets_[i].mergedScanId) == clusterBin) {
        unsigned int posInCluster = 0;
        BOOST_FOREACH (ScanId scannr, combineSets_[i].scans) {
          unsigned int fileIdx = fileList_.getFileIdx(scannr);
          //std::cerr << scannr << " " << i << " " << fileIdx << std::endl;
          SpectrumListPtr sl = msdVector[fileIdx / numMSFilePtrsPerBatch_]->run.spectrumListPtr;
          size_t result = getSpectrumIdxFromScannr(sl, hash_value(scannr));
          if (result >= sl->size()) {
            std::cerr << "  Warning: index " << result << " out of bounds: "
                  << fileList_.getFilePath(scannr) << ": " << scannr << std::endl;
          } else {
            scanIndices[fileIdx / numMSFilePtrsPerBatch_].push_back( MergeScanIndex(result, posInCluster, spectra.size()) );
            ++posInCluster;
          }
        }
        // initialize container for spectra to be merged for mergedScanId
        spectra.push_back(std::make_pair(std::vector<SpectrumPtr>(posInCluster), combineSets_[i].mergedScanId));
      }
    }

    std::cerr << "  Processing consensus spectra." << std::endl;
    for (unsigned int j = 0; j < numBatches_; ++j) {
      std::cerr << "  Batch " << j+1 << "/" << numBatches_ << std::endl;
      SpectrumListPtr sl = msdVector[j]->run.spectrumListPtr;
      std::sort(scanIndices[j].begin(), scanIndices[j].end(), lessIndex);
      BOOST_FOREACH (MergeScanIndex msi, scanIndices[j]) {
        SpectrumPtr s = sl->spectrum(msi.spectrumIndex, true);
        spectra[msi.mergeIdx].first[msi.posInCluster] = s;
      }
    }
    std::cerr << "  Merging spectra" << std::endl;
  #pragma omp parallel for schedule(dynamic, 100)
    for (int k = 0; k < static_cast<int>(spectra.size()); ++k) {
      if (spectra[k].first.size() > 0) {
        mergeSpectraSetMSCluster(spectra[k].first, spectra[k].second, mergedSpectra);
      }
    }
  }

  for (unsigned int i = 0; i < numBatches_; ++i) {
    std::string partSpecOutFN = getPartFN(spectrumOutFN_, "part" +
          boost::lexical_cast<std::string>(clusterBin) + "_" +
          boost::lexical_cast<std::string>(i) + "_tmpfile");
    if (remove(partSpecOutFN.c_str()) != 0) {
      std::cerr << "Warning: Can't remove " << partSpecOutFN << ": "
                << strerror(errno) << std::endl;
    }
  }
}

void MSFileMerger::writeClusterBins(unsigned int batchIdx,
    std::vector<MSDataPtr>& msDataPtrs, std::vector<SpectrumListSimplePtr>& spectrumLists) {
#pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < static_cast<int>(numClusterBins_); ++i) {
    if (spectrumLists[i]->spectra.size() == 0) {
      continue; // in case there are so few spectra that not all bins are filled
    }
    
    size_t idx = 0;
    BOOST_FOREACH (SpectrumPtr& s, spectrumLists[i]->spectra) {
      s->index = idx++;
    }
    
    MSDataMerger msdMerged(msDataPtrs);
    //MSData msdMerged;
    msdMerged.id = msdMerged.run.id = "merged_spectra";
    msdMerged.run.spectrumListPtr = spectrumLists[i];
    //msdMerged.dataProcessingPtrs = dataProcessingPtrs_;
    //msdMerged.instrumentConfigurationPtrs = instrumentConfigurationPtrs_;

    std::string partSpecOutFN = getPartFN(spectrumOutFN_, "part" +
          boost::lexical_cast<std::string>(i) + "_" +
          boost::lexical_cast<std::string>(batchIdx) + "_tmpfile");

    writeMSData(msdMerged, partSpecOutFN);
  }
}

} /* namespace maracluster */
