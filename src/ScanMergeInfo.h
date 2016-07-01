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
 
#ifndef SCANMERGEINFO_H
#define SCANMERGEINFO_H

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <climits>

#include <boost/foreach.hpp>

#include "ScanId.h"

struct ScanMergeInfo {
  ScanMergeInfo(ScanId _scannr, double _score, bool _isDecoy, 
                unsigned int _charge, std::string& _peptide) : 
      scannr(_scannr), score(_score), isDecoy(_isDecoy), 
      charge(_charge), peptide(_peptide), retentionTime(0.0) {}
	
	ScanMergeInfo(unsigned int _scannr, double _score, bool _isDecoy, 
                unsigned int _charge, std::string& _peptide) : 
      scannr(0,_scannr), score(_score), isDecoy(_isDecoy), 
      charge(_charge), peptide(_peptide), retentionTime(0.0) {}
      
	ScanMergeInfo(ScanId _scannr) : scannr(_scannr), score(1.0), 
	    isDecoy(false), charge(0), peptide(""), retentionTime(0.0) {}
  
	ScanMergeInfo() : scannr(0,0), score(1.0), isDecoy(false), 
	                  charge(0), peptide(""), retentionTime(0.0) {}
  
  ScanMergeInfo(const double _score, const std::string& _peptide, 
                const double _retentionTime) : 
      scannr(0,0), score(_score), isDecoy(false), charge(0), 
      peptide(_peptide), retentionTime(_retentionTime) {}
	
  ScanId scannr;
  double score;
  double retentionTime;
  bool isDecoy;
  unsigned int charge;
  std::string peptide;
  
  double weight;
  double precursorMZ;
  
  inline static bool lowerScannr(const ScanMergeInfo& a, 
	                               const ScanMergeInfo& b) { 
    return (a.scannr < b.scannr); 
  }
};

/* map: idx -> ScanMergeInfo */
typedef std::map<ScanId, ScanMergeInfo> ScanPeptideMap;

class ScanMergeInfoSet {
	public:
	  ScanMergeInfoSet() : isDecoy(false), minScannr_(UINT_MAX, UINT_MAX) {}
	  
	  inline size_t size() { return scans.size(); }
		void push_back(ScanMergeInfo smi) { 
		  scans.push_back(smi);
		  if (smi.scannr < minScannr_) minScannr_ = smi.scannr;
		}
		std::vector<ScanMergeInfo> scans;
		bool isDecoy;
		ScanId mergedScanId;
		double precursorMZ;
		std::string peptide;
		ScanId minScannr_;
		
		inline static bool lowerScore(const ScanMergeInfo& a, 
		                              const ScanMergeInfo& b) { 
      return (a.score < b.score); 
    }
		
		inline static bool lowerScannr(const ScanMergeInfoSet& a, 
		                              const ScanMergeInfoSet& b) { 
      return (a.minScannr_ < b.minScannr_); 
    }
    
		void sortByScore() {
		  std::sort(scans.begin(), scans.end(), lowerScore);
		}
};

std::ostream& operator<<(std::ostream& os, const ScanMergeInfoSet& sms);
/* map: peptide -> ScanMergeInfoSet */
typedef std::map<std::string, ScanMergeInfoSet> SpectraMergeMap;

#endif //SCANMERGEINFO_H
