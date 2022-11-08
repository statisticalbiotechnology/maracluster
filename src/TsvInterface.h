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
 
#ifndef MARACLUSTER_TSVINTERFACE_H_
#define MARACLUSTER_TSVINTERFACE_H_

#include <vector>
#include <cerrno>
#include <fstream>
#include <iostream>

#include <boost/foreach.hpp>

#include "MyException.h"

namespace maracluster {

// TODO: add exception handling
class TsvInterface {
 public:
  template <typename Type>
  static void write(const std::vector<Type>& vec, const std::string& outputFN, 
                    bool append) {
    std::ofstream outfile;
    if (append) {
      outfile.open(outputFN.c_str(), std::ios_base::app);
    } else {
      outfile.open(outputFN.c_str(), std::ios_base::out);
    }
    if (outfile.is_open()) {
      if (vec.size() > 0) {
        BOOST_FOREACH (const Type& t, vec) {
          outfile << t << '\n';
        }
      }
    } else {
      std::cerr << "(TsvInterface.h) could not open file for writing " << outputFN << std::endl;
    }  
  }
  
  template <typename Type>
  static void read(const std::string& inputFN, std::vector<Type>& vec) {
    if (!fileIsEmpty(inputFN)) {
      std::ifstream infile(inputFN.c_str(), std::ios_base::in);
      if (infile.is_open()) {
        while (infile.good()) {
          Type t;
          infile >> t;
          if (infile.good()) vec.push_back(t);
        }
      } else {
        std::cerr << "(TsvInterface.h) could not open file for reading " << inputFN << std::endl;
      }
    }
  }
  
  static bool fileIsEmpty(const std::string& fileName) {
    std::ifstream in(fileName.c_str(), std::ios::ate | std::ios::binary);
    if (in.is_open()) {
      return static_cast<long long>(in.tellg()) == 0ll;
    } else {
      return true;
    }
  }
};

// copied from percolator/DataSet.h
// using char pointers is much faster than istringstream
class TabReader {
 public:
  TabReader(const std::string& line) : f_(line.c_str()), err(0) {
    errno = 0;
  }
  
  void advance(const char* next) {
    if (*next != '\0') {
      f_ = next + 1; // eats up the tab
    } else {
      f_ = next; // prevents pointing over the null byte
    }
  }
  
  void skip() {
    const char* pch = strchr(f_, '\t');
    if (pch == NULL) {
      err = 1;
    } else {
      advance(pch);
    }
  }
  
  void skip(size_t numSkip) {
    for (size_t i = 0; i < numSkip; ++i) skip();
  }
  
  double readDouble() {
    char* next = NULL;
    errno = 0;
    double d = strtod(f_, &next);
    if (next == f_ || (*next != '\0' && !isspace(*next))
                   || ((d == HUGE_VAL || d == -HUGE_VAL) && errno == ERANGE)) {
      err = errno ? errno : 1;
    }
    advance(next);
    return d;
  }
  
  int readInt() {
    char* next = NULL;
    errno=0;
    long val = strtol(f_, &next, 10);
    if (next == f_ || (*next != '\0' && !isspace(*next))
                   || val < INT_MIN || val > INT_MAX) {
      err = errno ? errno : 1;
    }
    advance(next);
    return static_cast<int>(val);
  }
  
  std::string readString() {
    const char* pch = strchr(f_, '\t');
    if (pch == NULL) {
      err = 1;
      return std::string(f_);
    } else {
      std::string s(f_, static_cast<std::basic_string<char>::size_type>(pch - f_));
      advance(pch);
      return s;
    }
  }

  bool error() { return err != 0; }
 private:
  const char* f_;
  int err;
};


} /* namespace maracluster */

#endif /* MARACLUSTER_TSVINTERFACE_H_ */
