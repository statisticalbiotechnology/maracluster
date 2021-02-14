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
 
#ifndef MARACLUSTER_PVALBUFFER_H_
#define MARACLUSTER_PVALBUFFER_H_

#include <vector>
#include <string>
#include <cstdio>
#include <iterator>

#include <boost/foreach.hpp>

#include "BinaryInterface.h"
#include "MatrixLoader.h"
#include "PvalueTriplet.h"
#include "PvalueFilterAndSort.h"

namespace maracluster {

class PvalBuffer {
 public:
  // https://stackoverflow.com/questions/27593149/cannot-get-boost-foreach-to-work-with-my-custom-class
  template <typename Type>
  class PvalBufferIterator : public std::iterator<std::forward_iterator_tag, Type> {
   public:
    PvalBufferIterator(const PvalBufferIterator<Type>& pvalBufferIterator) {}
    
    PvalBufferIterator(const PvalBuffer& pvalBuffer, bool end = false) {
      if (end) {
        idx = -1;
      } else {
        if (pvalBuffer.isFile()) {
          matrixLoader_.initStream(pvalBuffer.getFileName());
        } else {
          matrixLoader_.initVector(pvalBuffer.getVector());
        }
      }
    }
    PvalBufferIterator& operator++() {
      idx++;
      if (idx >= pvalTriplets_.size()) {
        if (matrixLoader_.hasEdgesAvailable()) {
          pvalTriplets_.clear();
          matrixLoader_.nextNEdges(20e6, pvalTriplets_);
          idx -= 20e6;
        } else {
          idx = -1;
        }
      }
      return *this;
    }

    bool operator==(const PvalBufferIterator& other) const { return idx == other.idx; }
    bool operator!=(const PvalBufferIterator& other) const { return idx != other.idx; }
    Type& operator*() {
      if (idx == -1) {
        return nullValue_;
      } else {
        return pvalTriplets_.at(idx);
      }
    }
   private:
    MatrixLoader matrixLoader_;
    std::vector<PvalueTriplet> pvalTriplets_;
    int idx = 0;
    mutable PvalueTriplet nullValue_;
  };
  
  typedef PvalBufferIterator<PvalueTriplet> iterator_type;
  typedef PvalBufferIterator<const PvalueTriplet> const_iterator_type;
  
  iterator_type begin() { 
    return PvalBufferIterator<PvalueTriplet>(*this); 
  }
  iterator_type end() {
    return PvalBufferIterator<PvalueTriplet>(*this,true);
  }
  const_iterator_type begin() const {
    return PvalBufferIterator<const PvalueTriplet>(*this);
  }
  const_iterator_type end() const {
    return PvalBufferIterator<const PvalueTriplet>(*this ,true);
  }
    
  bool isFile() const {
    return !pvalFile.empty();
  }
  std::string getFileName() const {
    return pvalFile;
  }
  const std::vector<PvalueTriplet>& getVector() const {
    if (isFile()) {
      return pvalVec; //TODO: implement this    
    } else {
      return pvalVec;
    }
  }
  std::vector<PvalueTriplet>& getVector() {
    if (isFile()) {
      return pvalVec; //TODO: implement this    
    } else {
      return pvalVec;
    }
  }
  size_t size() const {
    if (isFile()) {
      return pvalCount;
    } else {
      return pvalVec.size();
    }
  }
  void clear() {
    if (isFile()) {
      remove(pvalFile.c_str());
    } else {
      clearVector();
    }
  }
  void clearVector() {
    std::vector<PvalueTriplet> empty;
    pvalVec.swap(empty);
  }
  void filterAndSort() {
    if (isFile()) {
      PvalueFilterAndSort::filterAndSort(pvalFile);
    } else {
      PvalueFilterAndSort::filterAndSort(pvalVec);
    }
  }
  void merge(PvalBuffer& other) {
    pvalCount += other.size();
    
    if (pvalCount > 40e6 || (other.isFile() && !isFile())) {
      bool append = false;
      pvalFile = std::tmpnam(nullptr);
      BinaryInterface::write<PvalueTriplet>(pvalVec, pvalFile, append);
      clearVector();
    }
    
    if (other.isFile()) {
      BinaryInterface::appendBinary(pvalFile, other.getFileName());
    } else {
      if (isFile()) { // append to current file
        bool append = true;
        BinaryInterface::write<PvalueTriplet>(other.getVector(), pvalFile, append);
      } else {
        pvalVec.insert(pvalVec.end(), other.getVector().begin(), other.getVector().end());
      }
    }
    other.clear();
  }
 private:
  std::vector<PvalueTriplet> pvalVec;
  std::string pvalFile;
  size_t pvalCount;
};

} /* namespace maracluster */

namespace boost {
  template<> struct range_mutable_iterator<maracluster::PvalBuffer> {
    typedef maracluster::PvalBuffer::iterator_type type;
  };

  template<> struct range_const_iterator<maracluster::PvalBuffer> {
    typedef maracluster::PvalBuffer::const_iterator_type type;
  };
}

#endif /* MARACLUSTER_PVALBUFFER_H_ */
