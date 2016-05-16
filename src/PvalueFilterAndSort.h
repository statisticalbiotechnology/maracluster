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
 
#ifndef PVALUE_FILTER_AND_SORT_H
#define PVALUE_FILTER_AND_SORT_H

#include <vector>
#include <map>
#include <set>
#include <string>

#include <iostream>
#include <sstream>
#include <fstream>
#include <cerrno>
#include <cstdio>
#include <ctime>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include <boost/filesystem.hpp>

#include "PvalueTriplet.h"
#include "BinaryInterface.h"

class PvalueFilterAndSort {
 public:
  static int maxPvalsPerFile_;
  
  static void filter(std::vector<PvalueTriplet>& buffer);
  static void filterAndSort(std::vector<PvalueTriplet>& buffer);
  static void filterAndSort(const std::string& pvalFN);
  static void filterAndSort(const std::vector<std::string>& pvalFNs, 
      const std::string& resultFN, bool tsvInput, bool removeUnidirected);
  
  static void filterAndSortSingleFile(const std::string& pvalFN, 
                                      bool removeUnidirected);
  
  static void externalMergeSort(const std::string& resultFN, int numFiles);
  
  static void convertBinaryPvalToTsv(std::string& binaryPvalFN, 
                                     std::string& tsvPvalFN);
  static bool unitTest();
  static bool singleFileUnitTest();
  
  inline static bool uniDirectionPval(const PvalueTriplet& a, 
                                      const PvalueTriplet& b) { 
    return (std::min(a.scannr1, a.scannr2) < std::min(b.scannr1, b.scannr2)) 
        || (std::min(a.scannr1, a.scannr2) == std::min(b.scannr1, b.scannr2) 
            && std::max(a.scannr1, a.scannr2) < std::max(b.scannr1, b.scannr2)) 
        || (std::min(a.scannr1, a.scannr2) == std::min(b.scannr1, b.scannr2) 
            && std::max(a.scannr1, a.scannr2) == std::max(b.scannr1, b.scannr2) 
            && a.pval > b.pval); 
  }
  
  inline static bool duplicatePval(const PvalueTriplet& a, 
                                   const PvalueTriplet& b) { 
    return (a.scannr1 < b.scannr1) || 
        (a.scannr1 == b.scannr1 && a.scannr2 < b.scannr2) || 
        (a.scannr1 == b.scannr1 && a.scannr2 == b.scannr2 && a.pval < b.pval); 
  }
  
  inline static bool scannr(const PvalueTriplet& a, const PvalueTriplet& b) { 
    return (a.scannr1 < b.scannr1) || 
        (a.scannr1 == b.scannr1 && a.scannr2 < b.scannr2); 
  }
                         
 private:
  static void reportProgress(time_t& startTime, clock_t& startClock);
  static int splitByHash(const std::vector<std::string>& pvalFNs, 
      const std::string& resultFN, bool tsvInput);
  
  static void readPvals(const std::string& pvalFN, 
                        std::vector<PvalueTriplet>& pvec);
  
  static void writePvalsTsv(std::ofstream& outfile, 
                            std::vector<PvalueTriplet>& pvec);
  
  static void removeDirectedDuplicates(std::vector<PvalueTriplet>& pvecIn, 
                                       std::vector<PvalueTriplet>& pvecOut);
  static void removeUndirectedDuplicates(std::vector<PvalueTriplet>& pvecIn, 
                                         std::vector<PvalueTriplet>& pvecOut);
  
  static long long estimateNumPvals(const std::vector<std::string>& pvalFNs, 
                                    bool tsvInput);
  static long long estimateNumPvals(const std::string& pvalFN, bool tsvInput);
  static long long getFileSize(const std::string& pvalFN);
  
  static void writeBufferToPartFiles(std::vector<PvalueTriplet>& buffer, 
                                     int numFiles, const std::string& resultFN);
};

#endif // PVALUE_FILTER_AND_SORT_H
