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
    const std::string& scanInfoFN, const std::string& scanDescFN,
    const std::string& resultBaseFN) {
  std::vector<PvalueTriplet> pvals;
  BOOST_FOREACH (const std::string& pvalTreeFN, pvalTreeFNs) {
    readPvalTree(pvalTreeFN, pvals);
  }
  std::sort(pvals.begin(), pvals.end());
  
  createScanDescriptionMap(scanInfoFN, scanDescFN, fileList);
  
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

void BatchSpectrumClusters::createScanDescriptionMap(
    const std::string& scanInfoFN, const std::string& scanDescFN,
    SpectrumFileList& fileList) {
  readScanNrs(scanInfoFN);
  readScanDescs(scanDescFN, fileList);
}

void BatchSpectrumClusters::readScanDescs(const std::string& scanDescFN,
    SpectrumFileList& fileList) {
  if (Globals::VERB > 1) {
    std::cerr << "Reading in scan descriptions." << std::endl;
  }
  if (!Globals::fileIsEmpty(scanDescFN)) {
    boost::iostreams::mapped_file mmap(scanDescFN, 
              boost::iostreams::mapped_file::readonly);
    const char* f = mmap.const_data();
    const char* l = f + mmap.size();
    
    errno = 0;
    char* next = NULL;
    const char* nextConst = NULL;
    unsigned int globalScannrTmp, localScannr;
    ScanId globalScannr;
    std::string filePath, peptide;
    double qval;
    
    bool hasScannrs = false;
    if (scanPeptideMap_.size() > 0) hasScannrs = true;
    
    while (errno == 0 && f && f<(l-12) ) {
      // we cannot trust this global scannr, calculate the correct one later
      globalScannrTmp = strtoul(f, &next, 0); f = next; 
      nextConst = strchr(f+1, '\t'); 
      filePath = std::string(f+1,nextConst); f = nextConst;
      localScannr = strtoul(f, &next, 0); f = next;
      nextConst = strchr(f+1, '\t'); 
      peptide = std::string(f+1,nextConst); f = nextConst;
      qval = strtod(f, &next); f = next;
      
      globalScannr = fileList.getScanId(filePath, localScannr);
      
      if (hasScannrs && scanPeptideMap_.find(globalScannr) != scanPeptideMap_.end()) {
        scanPeptideMap_[globalScannr] = 
            ScanMergeInfo(globalScannr, qval, false, 0, peptide);
      }
    }
  } else {
    std::cerr << "WARNING: Could not find scan desc file." << std::endl;
  }
}

void BatchSpectrumClusters::readScanNrs(const std::string& scanInfoFN) {
  if (Globals::fileExists(scanInfoFN)) {
    std::vector<ScanInfo> scanInfos;
    BinaryInterface::read<ScanInfo>(scanInfoFN, scanInfos);
    
    BOOST_FOREACH (const ScanInfo& si, scanInfos) {
      scanPeptideMap_[si.scanId] = ScanMergeInfo(si.scanId);
    }
  } else {
    std::cerr << "WARNING: Could not find scanInfo file." << std::endl;
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
      int negIntThreshold = 
          -1*static_cast<int>(clusterThresholds.at(thresholdIdx));
      std::string resultFN = resultBaseFN + "p" +
          boost::lexical_cast<std::string>(negIntThreshold) + ".tsv";
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
          << " " << clusters[pvalTriplet.scannr2].size() << " " << pvalTriplet.scannr2 << " " << pvalTriplet.pval << std::endl;   
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
    int negIntThreshold = 
        -1*static_cast<int>(clusterThresholds.at(thresholdIdx));
    std::string resultFN = resultBaseFN + "p" +
        boost::lexical_cast<std::string>(negIntThreshold) + ".tsv";
    writeClusters(clusters, fileList, resultFN);
    ++thresholdIdx;
  }
  
  if (Globals::VERB > 1) {
    std::cerr << "Finished writing clusterings." << std::endl;
  }
}

void BatchSpectrumClusters::writeClusters(
    std::map<ScanId, std::vector<ScanId> >& clusters,
    SpectrumFileList& fileList, const std::string& resultFN) {
  if (Globals::VERB > 2) {
    std::cerr << "Writing clustering to " << resultFN << std::endl;
  }
  std::ofstream resultStream(resultFN.c_str());
  std::map<ScanId, std::vector<ScanId> >::const_iterator it;
  std::set<ScanId> seenScannrs;
  for (it = clusters.begin(); it != clusters.end(); ++it) {
    if (it->second.size() > 0) {
      std::vector<ScanId>::const_iterator it2;
      for (it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
        ScanId globalScannr = *it2;
        seenScannrs.insert(globalScannr);
        
        unsigned int localScannr = fileList.getScannr(globalScannr);
        std::string filePath = fileList.getFilePath(globalScannr);
        std::string peptide = scanPeptideMap_[globalScannr].peptide;
        double qvalue = scanPeptideMap_[globalScannr].score;
        
        resultStream << filePath << '\t' << localScannr 
                     << '\t' << peptide << '\t' << qvalue << '\n';
      }
      resultStream << std::endl;
    }
  }
  writeSingletonClusters(seenScannrs, fileList, resultStream);
}
  
void BatchSpectrumClusters::writeSingletonClusters(
    std::set<ScanId>& seenScannrs, SpectrumFileList& fileList,
    std::ofstream& resultStream) {
  std::map<ScanId, ScanMergeInfo>::const_iterator spmIt;
  for (spmIt = scanPeptideMap_.begin(); spmIt != scanPeptideMap_.end(); ++spmIt) {
    if (seenScannrs.find(spmIt->first) == seenScannrs.end()) {
      ScanId globalScannr = spmIt->first;
      
      unsigned int localScannr = fileList.getScannr(globalScannr);
      std::string filePath = fileList.getFilePath(globalScannr);
      std::string peptide = spmIt->second.peptide;
      double qvalue = spmIt->second.score;
      resultStream << filePath << '\t' << localScannr 
                   << '\t' << peptide << '\t' << qvalue << '\n' << '\n';
    }
  }
}

bool BatchSpectrumClusters::scanDescReadUnitTest() {
  std::string scanInfoFN = "";
  std::string scanDescFN = "/media/storage/mergespec/data/103111-Yeast-2hr/percolator_no_tdc/scandesc/103111-Yeast-2hr.scannr_list.tsv";
  
  BatchSpectrumClusters clustering;
  SpectrumFileList fileList;
  clustering.scanPeptideMap_[ScanId(0,11)] = ScanMergeInfo(ScanId(0,11));
  clustering.scanPeptideMap_[ScanId(2,39214)] = ScanMergeInfo(ScanId(2,39214));
  clustering.createScanDescriptionMap(scanInfoFN, scanDescFN, fileList);
  if (clustering.scanPeptideMap_[ScanId(0,11)].peptide != "R.SIVPSGASTGVHEALEMR.D") {
    std::cerr << clustering.scanPeptideMap_[ScanId(0,11)].peptide << " != " 
              << "R.SIVPSGASTGVHEALEMR.D" << std::endl;
    return false;
  } else if (clustering.scanPeptideMap_[ScanId(2,39214)].peptide != "R.HSEFVAYPIQLLVTK.E") {
    std::cerr << clustering.scanPeptideMap_[ScanId(2,39214)].peptide << " != " 
              << "R.HSEFVAYPIQLLVTK.E" << std::endl;
    return false;
  } else {
    return true;
  } 
}

} /* namespace maracluster */
