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
 
#include "BatchSpectrumClusters.h"

namespace maracluster {

void BatchSpectrumClusters::printClusters(
    const std::vector<std::string>& pvalTreeFNs,
    const std::vector<double>& clusterThresholds, SpectrumFileList& fileList, 
    const std::string& scanInfoFN, const std::string& resultBaseFN) {
  std::vector<PvalueTriplet> pvals;
  BOOST_FOREACH (const std::string& pvalTreeFN, pvalTreeFNs) {
    readPvalTree(pvalTreeFN, pvals);
  }
  std::sort(pvals.begin(), pvals.end());
  
  if (Globals::fileExists(scanInfoFN)) {
    BinaryInterface::read<ScanInfo>(scanInfoFN, scanInfos_);
    std::sort(scanInfos_.begin(), scanInfos_.end());
  } else if (Globals::VERB > 1) {
    std::cerr << "WARNING: Could not find scanInfo file." << std::endl;
  }
  
  createClusterings(pvals, clusterThresholds, fileList, resultBaseFN);
}

void BatchSpectrumClusters::readPvalTree(const std::string& pvalTreeFN,
    std::vector<PvalueTriplet>& pvals) {
  if (Globals::VERB > 1) {
    std::cerr << "Reading in p-value tree." << std::endl;
  }
  if (!Globals::fileIsEmpty(pvalTreeFN)) {
    boost::iostreams::mapped_file mmap(pvalTreeFN, 
              boost::iostreams::mapped_file::readonly);
    const char* f = mmap.const_data();
    const char* l = f + mmap.size();
    
    errno = 0;
    char* next = NULL;
    PvalueTriplet tmp;
    while (errno == 0 && f && f <= (l-sizeof(tmp)) ) {
      tmp.readFromString(f, &next); f = next;
      pvals.push_back(tmp);
    }
  } else {
    std::cerr << "WARNING: Empty pvalue tree file " << pvalTreeFN << std::endl;
  }
}

void BatchSpectrumClusters::createClusterings(
    std::vector<PvalueTriplet>& pvals, 
    const std::vector<double>& clusterThresholds, SpectrumFileList& fileList,
    const std::string& resultBaseFN) {
  if (Globals::VERB > 1) {
    std::cerr << "Writing clusterings for " << clusterThresholds.size()
        << " thresholds." << std::endl;
  }
  unsigned int thresholdIdx = 0;
  std::map<ScanId, std::vector<ScanId> > clusters;
  std::map<ScanId, ScanId> mergeRoots;
  
  BOOST_FOREACH (PvalueTriplet& pvalTriplet, pvals) {
    while (pvalTriplet.pval > clusterThresholds.at(thresholdIdx)) {
      std::string resultFN = getClusterFN(resultBaseFN, clusterThresholds.at(thresholdIdx));
      writeClusters(clusters, fileList, resultFN);
      if (++thresholdIdx >= clusterThresholds.size()) break;
    }
    if (thresholdIdx >= clusterThresholds.size()) break;
    
    //std::cerr << pvalTriplet.scannr1 << " " << pvalTriplet.scannr2 << std::endl;
    if (mergeRoots.find(pvalTriplet.scannr1) == mergeRoots.end()) {
      clusters[pvalTriplet.scannr1].push_back(pvalTriplet.scannr1);
      mergeRoots[pvalTriplet.scannr1] = pvalTriplet.scannr1;
    } else {
      pvalTriplet.scannr1 = mergeRoots[pvalTriplet.scannr1];
    }
    
    if (mergeRoots.find(pvalTriplet.scannr2) == mergeRoots.end()) {
      clusters[pvalTriplet.scannr2].push_back(pvalTriplet.scannr2);
      mergeRoots[pvalTriplet.scannr2] = pvalTriplet.scannr2;
    } else {
      pvalTriplet.scannr2 = mergeRoots[pvalTriplet.scannr2];
    }
    
    if (pvalTriplet.scannr2 < pvalTriplet.scannr1) {
      std::swap(pvalTriplet.scannr1, pvalTriplet.scannr2);
    }
    
    if (clusters[pvalTriplet.scannr1].size() == 0) {
      std::cerr << "ERROR: Empty clusters: " << pvalTriplet.scannr1 
          << " " << clusters[pvalTriplet.scannr2].size() << " " 
          << pvalTriplet.scannr2 << " " << pvalTriplet.pval << std::endl;   
    }
    
    clusters[pvalTriplet.scannr1].insert(
        clusters[pvalTriplet.scannr1].end(),
        clusters[pvalTriplet.scannr2].begin(),
        clusters[pvalTriplet.scannr2].end());
    BOOST_FOREACH (ScanId scannr, clusters[pvalTriplet.scannr2]) {
      mergeRoots[scannr] = pvalTriplet.scannr1;
    }
    clusters.erase(pvalTriplet.scannr2);
  }
  
  while (thresholdIdx < clusterThresholds.size()) {
    std::string resultFN = getClusterFN(resultBaseFN, clusterThresholds.at(thresholdIdx));
    writeClusters(clusters, fileList, resultFN);
    ++thresholdIdx;
  }
  
  if (Globals::VERB > 1) {
    std::cerr << "Finished writing clusterings." << std::endl;
  }
}

std::string BatchSpectrumClusters::getClusterFN(
    const std::string resultBaseFN, double threshold) {
  int negIntThreshold = -1*static_cast<int>(threshold);
  std::string resultFN = resultBaseFN + "p" +
      boost::lexical_cast<std::string>(negIntThreshold) + ".tsv";
  return resultFN;
}

void BatchSpectrumClusters::writeClusters(
    std::map<ScanId, std::vector<ScanId> >& clusters,
    SpectrumFileList& fileList, const std::string& resultFN) {
  if (Globals::VERB > 2) {
    std::cerr << "Writing clusters to " << resultFN << std::endl;
  }
  std::ofstream resultStream(resultFN.c_str());
  std::map<ScanId, std::vector<ScanId> >::const_iterator it;
  std::set<ScanId> seenScannrs;
  std::vector<std::pair<size_t, size_t> > clusterSizeCounts(10);
  size_t clusterIdx = 1u;
  for (it = clusters.begin(); it != clusters.end(); ++it) {
    if (it->second.size() > 0) {
      size_t clusterSizeBin = 0u;
      size_t clusterSize = it->second.size();
      while (clusterSize >>= 1 && clusterSizeBin < 9) ++clusterSizeBin;
      ++clusterSizeCounts[clusterSizeBin].first;
      clusterSizeCounts[clusterSizeBin].second += it->second.size();
      
      std::vector<ScanId>::const_iterator it2;
      for (it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
        ScanId globalScannr = *it2;
        seenScannrs.insert(globalScannr);
        
        unsigned int localScannr = fileList.getScannr(globalScannr);
        std::string filePath = fileList.getFilePath(globalScannr);
        
        resultStream << filePath << '\t' << localScannr << '\t' << clusterIdx
                     << '\n';
      }
      clusterIdx++;
      resultStream << std::endl;
    }
  }
  size_t addedSingletons = writeSingletonClusters(seenScannrs, fileList, 
                                                  resultStream, clusterIdx);
  clusterSizeCounts[0].first += addedSingletons;
  clusterSizeCounts[0].second += addedSingletons;
  
  if (Globals::VERB > 2) {
    writeClusterSummary(clusterSizeCounts);
  }
}

void BatchSpectrumClusters::writeClusterSummary(
    std::vector<std::pair<size_t, size_t> >& clusterSizeCounts) {
  std::cerr << "clust_size\t#clusters\t#spectra" << std::endl;
  size_t totalClusters = 0u, totalSpectra = 0u;
  for (size_t i = 0u; i < 10u; ++i) {
    size_t minClusterSize = 1 << i;
    size_t maxClusterSize = (1 << (i+1)) - 1;
    if (i == 9u) {
      std::cerr << minClusterSize << "+";
    } else if (minClusterSize < maxClusterSize) {
      std::cerr << minClusterSize << "-" << maxClusterSize;
    } else {
      std::cerr << minClusterSize;
    }
    std::cerr << '\t' << clusterSizeCounts.at(i).first 
              << '\t' << clusterSizeCounts.at(i).second << std::endl;
    totalClusters += clusterSizeCounts.at(i).first;
    totalSpectra += clusterSizeCounts.at(i).second;
  }
  std::cerr << "total\t" << totalClusters << '\t' << totalSpectra << std::endl;
  std::cerr << std::endl;
}
  
size_t BatchSpectrumClusters::writeSingletonClusters(
    std::set<ScanId>& seenScannrs, SpectrumFileList& fileList,
    std::ofstream& resultStream, size_t clusterIdx) {
  size_t addedSingletonClusters = 0u;
  std::vector<ScanInfo>::const_iterator spmIt;
  for (spmIt = scanInfos_.begin(); spmIt != scanInfos_.end(); ++spmIt) {
    if (seenScannrs.find(spmIt->scanId) == seenScannrs.end()) {
      ScanId globalScannr = spmIt->scanId;
      
      addedSingletonClusters += 1;
      
      unsigned int localScannr = fileList.getScannr(globalScannr);
      std::string filePath = fileList.getFilePath(globalScannr);
      resultStream << filePath << '\t' << localScannr << '\t' << clusterIdx++
                   << '\n' << '\n';
    }
  }
  return addedSingletonClusters;
}

} /* namespace maracluster */
