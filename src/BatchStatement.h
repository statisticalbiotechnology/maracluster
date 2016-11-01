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
 
#ifndef BATCH_STATEMENT_H
#define BATCH_STATEMENT_H

#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

#include "MyException.h"
#include "Globals.h"

enum BatchStatementType {
  kSpectra, kPvalueVectors, kPvalueVectorsHead, kPvalueVectorsTail,
  kPvalues, kMassRange, kGetPeakCounts, 
  kGetArchivedPeakCounts, kInsertPeakCountsArchive, kInsertPeakCounts,
  kGetFilePeakCounts, kFiles
};

class BatchStatement {
 public:
  BatchStatement(BatchStatementType t) : type(t) { }
  BatchStatement(BatchStatementType t, size_t s) : type(t) {
    strings.resize(s);
  }
  
  void execute() { }
  void execute(std::vector<BatchStatement>& result) { }
  void execute(std::ostream& outStream) {
    switch (type) {
      case kPvalueVectors:
      case kPvalueVectorsHead:
      case kPvalueVectorsTail:
      #pragma omp critical (write_pval_vec)
        {
          write(outStream);
        }
        break;
      default:
        break;
    }
  }
  
  void write(std::ostream& outStream) {
    BOOST_FOREACH (std::string& s, strings) {
      outStream << s << "\t";
    }
    outStream << "\n";
  }
                             
  void read(std::istream& inStream) {
    std::string col;
    while (std::getline(inStream, col, '\t')) {
      strings.push_back(col);
    }
  }
  
  inline void bindBuffer(std::stringstream *buf) {
    buffer = buf;
  }
  
  inline void bindInt32(const int pos, const int i) {
    strings[pos] = "i:" + boost::lexical_cast<std::string>(i);
  }
  
  inline void bindFloat(const int pos, const float f) {
    strings[pos] = "f:" + boost::lexical_cast<std::string>(f);
  }
  
  inline void bindDouble(const int pos, const double d) {
    strings[pos] = "d:" + boost::lexical_cast<std::string>(d);
  }
  
  inline void bindString(const int pos, const std::string& s) {
    strings[pos] = "s:" + s;
  }
  
  inline void bindBytes(const int pos, const std::string& b) {
    strings[pos] = "b:" + b;
  }
  
  inline int getInt32(const int pos) {
    if (strings[pos].size() > 0) {
      int i = 0;
      if (strings[pos].substr(0,1) != "i")
        std::cerr << "Warning: wrong type" << std::endl;
      i = atoi(strings[pos].substr(2).c_str());
      return i;
    } else {
      return 0;
    }
  }
  
  inline float getFloat(const int pos) {
    if (strings[pos].size() > 0) {
      float f = 0.0f;
      if (strings[pos].substr(0,1) != "f")
        std::cerr << "Warning: wrong type" << std::endl;
      f = static_cast<float>(atof(strings[pos].substr(2).c_str()));
      return f;
    } else {
      return 0.0f;
    }
  }
  
  inline double getDouble(const int pos) {
    if (strings[pos].size() > 0) {
      double d = 0.0;
      if (strings[pos].substr(0,1) != "d")
        std::cerr << "Warning: wrong type" << std::endl;
      d = atof(strings[pos].substr(2).c_str());
      return d;
    } else {
      return 0.0;
    }
  }
  
  inline std::string getString(const int pos) {
    if (strings[pos].size() > 0) {
      std::string s = "";
      if (strings[pos].substr(0,1) != "s")
        std::cerr << "Warning: wrong type" << std::endl;
      s = strings[pos].substr(2);
      return s;
    } else {
      return "";
    }
  }
  
  inline void getBytes(const int pos, std::string& b) {
    if (strings[pos].size() > 0) {
      if (strings[pos].substr(0,1) != "b")
        std::cerr << "Warning: wrong type" << std::endl;
      b = strings[pos].substr(2);
    }
  }
 private: 
  BatchStatementType type;
  std::vector<std::string> strings;
  std::stringstream *buffer;
};

#endif // BATCH_STATEMENT_H
