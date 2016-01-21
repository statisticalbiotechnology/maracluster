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
 
#include "PvalueFilterAndSort.h"

int PvalueFilterAndSort::maxPvalsPerFile_ = 5000000;

void PvalueFilterAndSort::filterAndSort(const std::string& pvalFN) {
  bool tsvInput = false;
  bool removeUnidirected = false;
  std::vector<std::string> pvalFNs;
  pvalFNs.push_back(pvalFN);
  filterAndSort(pvalFNs, pvalFN, tsvInput, removeUnidirected);
}

void PvalueFilterAndSort::filterAndSort(const std::vector<std::string>& pvalFNs, 
    const std::string& resultFN, bool tsvInput, bool removeUnidirected) {
  int numFiles = splitByHash(pvalFNs, resultFN, tsvInput);

#pragma omp parallel for schedule(dynamic, 1)  
  for (int bin = 0; bin < numFiles; ++bin) {
    if (bin % 10 == 0) {
      std::cerr << "Sorting and filtering bin " << bin+1 << "/" << numFiles << " (" << (bin+1)*100/numFiles << "%)" << std::endl;
    }
    std::string partFileFN = resultFN + "." + boost::lexical_cast<std::string>(bin);
    
    // ignore files that do not exist
    if (!boost::filesystem::exists(partFileFN)) { 
        std::cerr << "Ignoring missing result file \"" << 
                partFileFN << "\" does not exist." << std::endl;
        continue;
    }
    
    filterAndSortSingleFile(partFileFN, removeUnidirected);
  }
  
  externalMergeSort(resultFN, numFiles);  
}

int PvalueFilterAndSort::splitByHash(const std::vector<std::string>& pvalFNs,
    const std::string& resultFN, bool tsvInput) {
  long long numPvals = estimateNumPvals(pvalFNs, tsvInput);
  int numFiles = numPvals/maxPvalsPerFile_ + 1;
  
  if (numPvals == 0) return 0;
  
  std::cerr << "Estimated " << numPvals << " p-values in this file" << std::endl;
  std::cerr << "Writing " << numFiles << " part files" << std::endl;
  
  std::vector<PvalueTriplet> buffer;
  buffer.reserve(maxPvalsPerFile_);
  
  long long i = 0;
  BOOST_FOREACH (const std::string& pvalFN, pvalFNs) {
    if (estimateNumPvals(pvalFN, tsvInput) == 0) continue;
    boost::iostreams::mapped_file mmap(pvalFN, 
            boost::iostreams::mapped_file::readonly);
    const char* f = mmap.const_data();
    const char* l = f + mmap.size();
    
    {
      errno = 0;
      char* next = NULL;
      PvalueTriplet tmp;
      while (errno == 0 && f && f <= (l-sizeof(tmp)) ) {
        if (tsvInput) {
          tmp.readFromString(f, &next); f = next;
        } else {
          memcpy(&tmp, f, sizeof(tmp));
          f += sizeof(tmp);
        }
        buffer.push_back(tmp);
        //std::cerr << tmp.scannr1 << " " << tmp.scannr2 << " " << tmp.pval << std::endl;
        if (++i % maxPvalsPerFile_ == 0) {
          std::cerr << "Hashing p-value " << i << " (" << i*100/numPvals << "%)" << std::endl;
          writeBufferToPartFiles(buffer, numFiles, resultFN);
          //std::cerr << "Written p-value " << i << "/" << numPvals << std::endl;
          buffer.clear();
          buffer.reserve(maxPvalsPerFile_);
        }
      }
    }
  }

  if (buffer.size() > 0) writeBufferToPartFiles(buffer, numFiles, resultFN);  
  return numFiles;
}

void PvalueFilterAndSort::filterAndSortSingleFile(const std::string& partFileFN,
    bool removeUnidirected) {
  std::vector<PvalueTriplet> buffer, filteredBuffer;
  buffer.reserve(maxPvalsPerFile_);
  
  readPvals(partFileFN, buffer);
  
  size_t numPvals = buffer.size();
  
  //std::cerr << "Read " << numPvals << " p-values" << std::endl;
  
  removeDirectedDuplicates(buffer, filteredBuffer);
  
  if (removeUnidirected) {
    removeUndirectedDuplicates(filteredBuffer, buffer);
  } else {
    buffer.swap(filteredBuffer);
  }
  
  //std::cerr << "Removed " << numPvals - buffer.size() << " duplicate entries" << std::endl;
  
  std::sort(buffer.begin(), buffer.end(), lowerPval);
  
  //std::string sortedPvalFN = partFileFN + ".sorted";
  std::string sortedPvalFN = partFileFN;
  
  bool append = false;
  BinaryInterface::write<PvalueTriplet>(buffer, sortedPvalFN, append);
}

void PvalueFilterAndSort::externalMergeSort(const std::string& resultFN, int numFiles) {
  //size_t pageSize = boost::mapped_region::get_page_size();
  
  int numOpenFiles = numFiles;  
  std::vector<PvalueTriplet> pvecBuffer;
  std::vector<long long> offsets(numFiles);
  int maxPvalSort = maxPvalsPerFile_;
  bool first = true;
  long long numPvals = 0, numWrittenPvals = 0;
  while (numOpenFiles > 0) {
    float maxMinPval = 0.0;
    
    for (int bin = 0; bin < numFiles; ++bin) {
      std::string partFileFN = resultFN + "." + boost::lexical_cast<std::string>(bin);
      
      // ignore files that do not exist
      if (!boost::filesystem::exists(partFileFN)) { 
        std::cerr << "Ignoring missing result file \"" << 
            partFileFN << "\" does not exist." << std::endl;
        continue;
      }
      
      boost::iostreams::mapped_file mmap(partFileFN, 
            boost::iostreams::mapped_file::readonly);
      
      const char* f = mmap.const_data();
      const char* l = f + mmap.size();
      f += offsets[bin];
      
      if (first) numPvals += (l - f)/sizeof(PvalueTriplet);
      
      if (f && f <= (l-sizeof(PvalueTriplet))) {
        errno = 0;
        PvalueTriplet tmp;
        int numPvalsAdded = 0;
        while (errno == 0 && f && f<=(l-sizeof(tmp)) && 
               ++numPvalsAdded <= maxPvalSort/numFiles) {
          memcpy(&tmp, f, sizeof(tmp));
          f += sizeof(tmp);
          offsets[bin] += sizeof(tmp);
          pvecBuffer.push_back(tmp);
        }
        
        if (pvecBuffer.back().pval < maxMinPval) {
          maxMinPval = pvecBuffer.back().pval;
        }
      }
      
      if (!(f && f <= (l-sizeof(PvalueTriplet)))) --numOpenFiles;
    }
    
    if (first) first = false;
    
    std::sort(pvecBuffer.begin(), pvecBuffer.end(), lowerPval);
    
    std::vector<PvalueTriplet> pvec;
    int numPvalsToErase = 0;
    BOOST_FOREACH (PvalueTriplet p, pvecBuffer) {
      if (numOpenFiles == 0 || p.pval <= maxMinPval) {
        pvec.push_back(p);
        ++numPvalsToErase;
      } else {
        break;
      }
    }
    pvecBuffer.erase(pvecBuffer.begin(), pvecBuffer.begin() + numPvalsToErase);
    
    numWrittenPvals += pvec.size();
    std::cerr << "Writing p-value " << numWrittenPvals << "/" << numPvals << " (" << numWrittenPvals*100/numPvals << "%)"<< std::endl;
    
    bool append = true;
    BinaryInterface::write<PvalueTriplet>(pvec, resultFN, append);
  }
  
  for (int bin = 0; bin < numFiles; ++bin) {
    std::string partFileFN = resultFN + "." + boost::lexical_cast<std::string>(bin);
    remove(partFileFN.c_str());
  }
}

void PvalueFilterAndSort::writeBufferToPartFiles(std::vector<PvalueTriplet>& buffer, int numFiles, const std::string& resultFN) {
  std::vector< std::vector<PvalueTriplet> > hashBins(numFiles);
  BOOST_FOREACH (PvalueTriplet t, buffer) {
    // TODO: this function might not work correctly
    int hashBin = ((t.scannr1.scannr % numFiles) + (t.scannr2.scannr % numFiles)) % numFiles;
    hashBins[hashBin].push_back(t);
  }
  
  //std::cerr << "Start writing p-values" << std::endl;
  
  int j = 0;
  typedef std::vector<PvalueTriplet> PvalTripletVector;
  BOOST_FOREACH (PvalTripletVector pvec, hashBins) {
    if (pvec.size() > 0) {
      std::string partFileFN = resultFN + "." + boost::lexical_cast<std::string>(j);
      bool append = true;
      BinaryInterface::write<PvalueTriplet>(pvec, partFileFN, append);
    }
    ++j;
  }
}

long long PvalueFilterAndSort::estimateNumPvals(const std::vector<std::string>& pvalFNs, bool tsvInput) {
  long long acc = 0;
  BOOST_FOREACH (const std::string& pvalFN, pvalFNs) {
    acc += estimateNumPvals(pvalFN, tsvInput);
  }
  return acc;
}

long long PvalueFilterAndSort::estimateNumPvals(const std::string& pvalFN, bool tsvInput) {
  long long fileSize = getFileSize(pvalFN);
  if (tsvInput) {
    std::ifstream pvalStream(pvalFN.c_str());
    std::string line;
    
    int sampleSize = 1000;
    int samplesRead = 0;
    double bytes = 0;
    if (pvalStream.is_open()) {
      while (getline(pvalStream, line)) {
        bytes += line.size();
        if (++samplesRead >= sampleSize) break;
      }
    }
    
    std::cerr << "Est. bytes per row: " << bytes / sampleSize << std::endl;
    std::cerr << "Est. p-values: " << fileSize / (bytes / sampleSize) << std::endl;
    return static_cast<long long>(fileSize / (bytes / sampleSize));
  } else {
    return static_cast<long long>(fileSize / sizeof(PvalueTriplet));
  }
}

long long PvalueFilterAndSort::getFileSize(const std::string& pvalFN) {
  std::ifstream in(pvalFN.c_str(), std::ios::ate | std::ios::binary);
  if (in.is_open()) {
    return static_cast<long long>(in.tellg());
  } else {
    std::cerr << "WARNING: could not read any p-values from "
        << pvalFN << "." << std::endl;
    return 0LL;
  }
}

// keep only the lowest p-value per directed pair
void PvalueFilterAndSort::removeDirectedDuplicates(
    std::vector<PvalueTriplet>& pvecIn, 
    std::vector<PvalueTriplet>& pvecOut) {
  std::sort(pvecIn.begin(), pvecIn.end(), duplicatePval);
  
  pvecOut.clear();
  pvecOut.reserve(pvecIn.size());
  
  ScanId lastScannr1, lastScannr2;
  for (std::vector<PvalueTriplet>::iterator it = pvecIn.begin(); it != pvecIn.end(); ++it) {
    if (!(it->scannr1 == lastScannr1 && it->scannr2 == lastScannr2)) {
      pvecOut.push_back(*it);
      lastScannr1 = it->scannr1;
      lastScannr2 = it->scannr2;
    }
  }  
}

// keep only the highest p-value per undirected pair
void PvalueFilterAndSort::removeUndirectedDuplicates(
    std::vector<PvalueTriplet>& pvecIn, 
    std::vector<PvalueTriplet>& pvecOut) {
  std::sort(pvecIn.begin(), pvecIn.end(), uniDirectionPval);
  
  pvecOut.clear();
  pvecOut.reserve(pvecIn.size());
  
  ScanId lastScannr1, lastScannr2;
  pvecOut.reserve(pvecIn.size());
  for (std::vector<PvalueTriplet>::iterator it = pvecIn.begin(); it != pvecIn.end(); ++it) {
    if (!(it->scannr1 == lastScannr2 && it->scannr2 == lastScannr1)) {
      pvecOut.push_back(*it);
      lastScannr1 = it->scannr1;
      lastScannr2 = it->scannr2;
    }
  }  
}

// TODO: check if file exists before reading
void PvalueFilterAndSort::readPvals(const std::string& pvalFN, std::vector<PvalueTriplet>& pvec) {
#pragma omp critical (read_pval_parts)
  {
    BinaryInterface::read<PvalueTriplet>(pvalFN, pvec);
  }
}

void PvalueFilterAndSort::convertBinaryPvalToTsv(std::string& binaryPvalFN, 
    std::string& tsvPvalFN) {
  std::vector<PvalueTriplet> pvec;
  readPvals(binaryPvalFN, pvec);
  std::ofstream outfile(tsvPvalFN.c_str(), std::ios::out);
  writePvalsTsv(outfile, pvec);
}

void PvalueFilterAndSort::writePvalsTsv(std::ofstream& outfile, std::vector<PvalueTriplet>& pvec) {
  if (outfile.is_open()) {
    BOOST_FOREACH (PvalueTriplet t, pvec) {
      outfile << t << "\n";
    }
  }
}

bool PvalueFilterAndSort::unitTest() {
  std::string pvalFN = "/media/storage/mergespec/data/batchtest/1300.ms2.pvalues_subset.tsv";
  //std::string pvalFN = "/media/storage/mergespec/data/batchcluster/Linfeng/Linfeng.pval_matrix_triplets.tsv";
  std::string resultFN = pvalFN + ".sorted.tsv";
 
  bool removeUnidirected = true;
  bool tsvInput = true;
  std::vector<std::string> pvalFNs;
  pvalFNs.push_back(pvalFN);
  filterAndSort(pvalFNs, resultFN, tsvInput, removeUnidirected);
  //externalMergeSort(pvalFN, 4);
  /*
  std::string binaryPvalFN = pvalFN + ".2";
  std::string tsvPvalFN = pvalFN + ".2.tsv";
  convertBinaryPvalToTsv(binaryPvalFN, tsvPvalFN);
  */
  return true;
}

bool PvalueFilterAndSort::singleFileUnitTest() {
  std::string pvalPartFN = "/media/storage/mergespec/data/batchtest/1300.ms2.pvalues_subset.tsv.3";
  bool removeUnidirected = true;
  filterAndSortSingleFile(pvalPartFN, removeUnidirected);
  
  return true;
}
