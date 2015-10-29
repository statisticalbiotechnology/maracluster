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
 
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <vector>

#include <boost/foreach.hpp>

#include "ScanMergeInfo.h"
#include "SpectrumFileList.h"

unsigned int extractScannr(std::string id) {
  std::vector< std::string > elems;
  std::stringstream ss(id);
  std::string item;
  while (std::getline(ss, item, '_')) {
    elems.push_back(item);
  }
  return static_cast<unsigned int>(std::strtoul(elems[elems.size() - 3].c_str(), NULL, 0));
}

std::string extractCharge(std::string id) {
  std::vector< std::string > elems;
  std::stringstream ss(id);
  std::string item;
  while (std::getline(ss, item, '_')) {
    elems.push_back(item);
  }
  return elems[elems.size() - 2];
}

int extractIntCharge(std::string id) {
  std::vector< std::string > elems;
  std::stringstream ss(id);
  std::string item;
  while (std::getline(ss, item, '_')) {
    elems.push_back(item);
  }
  return atoi(elems[elems.size() - 2].c_str());
}

/* Creates a map: ("peptide" + "_" + "charge") => ScanMergeInfo(scanNr, PEP, 
     isDecoy, charge, peptide) of all PSMs under qvalue_thresh */
void parsePercOutfile(std::string percOutFN, SpectraMergeMap& peptideScanMap,
                      bool isDecoy, double qvalue_thresh = 0.01) {
  std::ifstream percOut(percOutFN.c_str());
  std::string line, tmp, id, peptide;
  double qvalue, pep;
  if (percOut.is_open()) {
    getline(percOut, line); /* remove header line */
    while (getline(percOut, line)) {
      std::istringstream iss(line);
      iss >> id >> tmp >> qvalue >> pep >> peptide;
      std::string key = peptide + "_" + extractCharge(id);
      if (qvalue > qvalue_thresh) break;
      peptideScanMap[key].push_back(ScanMergeInfo(extractScannr(id), pep,
                                     isDecoy, extractIntCharge(id), peptide));
      /* std::cerr << "Inserting peptide " << peptide << " with scannr " <<
                       scannr << std::endl; */
    }
  } else {
    throw std::runtime_error("Could not open percolator output file");
  }
  percOut.close();
}

/* Creates a map: scannr => ScanMergeInfo(scanNr, qvalue, isDecoy, charge, 
     peptide) of all PSMs under qvalue_thresh */
void parsePercOutfile(std::string percOutFN, ScanPeptideMap& scanPeptideMap,
                      bool isDecoy = false) {
  std::ifstream percOut(percOutFN.c_str());
  std::string line, tmp, id, peptide;
  double qvalue, pep;
  if (percOut.is_open()) {
    getline(percOut, line); /* remove header line */
    while (getline(percOut, line)) {
      std::istringstream iss(line);
      iss >> id >> tmp >> qvalue >> pep >> peptide;
      unsigned int scannr = extractScannr(id);
      int charge = extractIntCharge(id);
      // TODO: Fix this!
      ScanId scanId(0,scannr);
      if (qvalue < scanPeptideMap[scanId].score) {
        scanPeptideMap[scanId] = 
            ScanMergeInfo(scannr, qvalue, isDecoy, charge, peptide);
      }
    }
  } else {
    throw std::runtime_error("Could not open percolator output file");
  }
  percOut.close();
}

/* Filters out the duplicates from the map produced by parsePercOutfile and 
   creates a vector of ScanMergeInfoSets */
void reducePercOutfile(SpectraMergeMap& peptideScanMap, 
    const std::string spectrumFileIn,
    SpectrumFileList& fileList, std::vector<ScanMergeInfoSet>& combineSets, 
    const std::string scanWeightsFN) {
  int count = 0;
  int decoyCount = 0;
  std::set< std::string > uniquePeptides;
  
  // decide if we write the merger list to a file or not at all
  std::ofstream scanWeightsStream;
  bool writeMergeInfo = false;
  if (scanWeightsFN.size() > 0) {
    scanWeightsStream.open(scanWeightsFN.c_str());
    writeMergeInfo = true;
  }
  
  typedef SpectraMergeMap::value_type PeptideScanPair;
  
  /* this ensures the right order of the files in the fileList object */
  ScanId tmp = fileList.getScanId(spectrumFileIn, 0);
  std::string mergedFilepath = "merged";
  std::string mergedDecoyFilepath = "merged_decoy";
  tmp = fileList.getScanId(mergedFilepath, 0);
  tmp = fileList.getScanId(mergedDecoyFilepath, 0);
  
  BOOST_FOREACH (PeptideScanPair& peptideScanPair, peptideScanMap) {
  	ScanMergeInfoSet& scanMergeInfoSet = peptideScanPair.second;
    if (scanMergeInfoSet.scans.size() >= 2) {
      bool first = true;
      double norm = 0;
      int decoyMult = 1;
      BOOST_FOREACH (ScanMergeInfo& scanMergeInfo, scanMergeInfoSet.scans) {
        if (first) {
          norm = scanMergeInfo.score;
          scanMergeInfo.weight = 1;
          first = false;
          std::string filepath = mergedFilepath;
          if (scanMergeInfo.isDecoy) {
		      	scanMergeInfoSet.isDecoy = true;
		      	filepath = mergedDecoyFilepath;
		      	++decoyCount;
		     	} else {
		     	  uniquePeptides.insert(peptideScanPair.first.substr(0,
		     	                            peptideScanPair.first.size() - 2));
		     	}
          scanMergeInfoSet.mergedScanId = 
              fileList.getScanId(filepath, ++count);
          scanMergeInfoSet.peptide = peptideScanPair.first;
        } else {
          scanMergeInfo.weight = norm/scanMergeInfo.score;
        }
      }
      if (writeMergeInfo) {
        scanWeightsStream << scanMergeInfoSet;
      }
      combineSets.push_back(scanMergeInfoSet);
    }
  }
  /* std::cout << "Found " << combineSets.size() << " peptide matches with " <<
     "multiple scan matches: " << combineSets.size() - decoyCount << 
     " target merges (" << uniquePeptides.size() << " unique peptides) and " <<
     decoyCount << " decoy merges" << std::endl; */
}

void reducePercInfile(std::string percInFN, std::string percInReducedFN, 
                      SpectraMergeMap& peptideScanMap) {
  std::ifstream percIn(percInFN.c_str());
  std::ofstream percInReduced(percInReducedFN.c_str());
  std::string line, tmp, id, peptide;
  int label;
  double precursorMZ;
  std::vector<double>::iterator it;
  std::vector<double> decoyPrecursorMasses;
  unsigned int scannr;
  if (percIn.is_open()) {
    getline(percIn, line); /* read header line */
    percInReduced << line << std::endl; /* print header line to file */
    
    std::istringstream issh(line);
    unsigned int index = 0, mindex, pindex;
    while (issh.good()) {
      issh >> tmp;
      if (tmp == "Mass") mindex = index;
      else if (tmp == "Peptide") pindex = index;
      index++;
    }
    if (mindex == 0 || pindex == 0) 
      throw std::runtime_error("Could not find Mass and Peptide headers");
    while (getline(percIn, line)) {
      std::istringstream iss(line);
      iss >> id >> label >> scannr;
      if (id == "DefaultDirection") continue;
      for (unsigned int i = 0; i < mindex - 3; ++i) iss >> tmp;
      iss >> precursorMZ;
      for (unsigned int i = 0; i < pindex - mindex - 1; ++i) iss >> tmp;
      iss >> peptide;
      std::string key = peptide + "_" + extractCharge(id);
      if (peptideScanMap[key].scans.size() >= 2) {
        peptideScanMap[key].precursorMZ = precursorMZ;
        /*std::cerr << "Adding precursorMZ " << precursorMZ << " to key " <<
                       key << std::endl; */
      	decoyPrecursorMasses.push_back(precursorMZ);
      	continue;
     	} else if (label < 0 && 
          (it = find(decoyPrecursorMasses.begin(), decoyPrecursorMasses.end(),
                     precursorMZ)) != decoyPrecursorMasses.end()) {
     		decoyPrecursorMasses.erase(it);
     		continue;
     	}
      percInReduced << line << std::endl;
      /*std::cerr << "Inserting peptide " << peptide << " with scannr " <<
        scannr << " and mass " << precursorMZ << std::endl; */
    }
  } else {
    throw std::runtime_error("Could not open percolator input file");
  }
  percIn.close();
  percInReduced.close();
}
