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
 
#ifndef MARACLUSTER_PVALUECALCULATOR_H_
#define MARACLUSTER_PVALUECALCULATOR_H_

/**
This file provides functions for calculating pvalues using dynamic programming
Adapted from https://code.google.com/p/in-silico-mass-fingerprinting/
*/
#include <vector>
#include <map>
#include <string>

#include <cmath>
#include <algorithm>
#include <numeric>

#include <iostream>
#include <sstream>

#include <cassert>
#include <cstdlib>
#include <stdexcept>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace maracluster {

class PvalueCalculator {
 public:
  static unsigned int probDiscretizationLevels_;
  static const unsigned int kPolyfitDegree = 5u;
  static const double kMinProb, kMaxProb;
  static const unsigned int kMaxScoringPeaks = 40u;
  static unsigned int kMinScoringPeaks;
  static const bool kVariableScoringPeaks;
  
  PvalueCalculator() : maxScore_(0u) {}
  inline unsigned int getNumScoringPeaks() const { return peakBins_.size(); }
  
  // prevents ODR usage of kMaxScoringPeaks: http://en.cppreference.com/w/cpp/language/definition#ODR-use
  static unsigned int getMaxScoringPeaksConstant() {
    return kMaxScoringPeaks;
  }
  
  static unsigned int getMaxScoringPeaks(double mass) { 
    unsigned int maxScoringPeaks = getMaxScoringPeaksConstant();
    if (kVariableScoringPeaks) {
      return std::min(maxScoringPeaks, static_cast<unsigned int>(mass / 50.0));
    } else {
      return maxScoringPeaks;
    }
  }
  
  static unsigned int getMinScoringPeaks(double mass) { 
    if (kVariableScoringPeaks) {
      return std::min(kMinScoringPeaks, static_cast<unsigned int>(mass / 50.0));
    } else {
      return kMinScoringPeaks;
    }
  }
  
  inline void setPeakBins(std::vector<unsigned int>& peakBins) {
    peakBins_.swap(peakBins);
  }
  void initPolyfit(std::vector<unsigned int>& peakBins, 
                   std::vector<unsigned int>& peakScores, 
                   std::vector<double>& polyfit);
  
  void computePvalVector(const std::vector<double>& peakDist, 
                         std::vector<double>& sumProb);
  double computePval(const std::vector<unsigned int>& queryPeakBins, 
      bool smoothing, std::vector<double>& sumProb);
  
  void computePvalVectorPolyfit(const std::vector<double>& peakDist);
  double computePvalPolyfit(const std::vector<unsigned int>& queryPeakBins);
  
  void copyPolyfit(short* peakBins, short* peakScores, double* polyfit);
  void serialize(std::string& polyfitString, std::string& peakScorePairsString);
  void deserialize(std::string& polyfitString, std::string& peakScorePairsString);
  
  inline std::vector<unsigned int> getPeakBins() const { return peakBins_; }
  inline std::vector<unsigned int>& getPeakBinsRef() { return peakBins_; }
  inline std::vector<unsigned int> getPeakScores() const { return peakScores_; }
  inline std::vector<double> getPolyfit() const { return polyfit_; }
  
  static bool pvalUnitTest();
  static bool pvalPolyfitUnitTest();
  static bool pvalUniformUnitTest();
  static bool binaryPeakMatchUnitTest();
  
  // needed for smoothing and unit tests
  inline static void setSeed(unsigned long s) { seed_ = s; }
  static unsigned long lcg_rand();
  static double lcg_rand_unif();
  
 private:
  unsigned int maxScore_;
  std::vector<unsigned int> peakBins_;
  std::vector<unsigned int> peakScores_;
  std::vector<double> polyfit_;
  
  void initFromPeakBins(const std::vector<double>& peakDist,
                        std::vector<double>& peakProbs);
  void binaryMatchPeakBins(const std::vector<unsigned int>& queryPeakBins, std::vector<bool>& d);
  double polyval(double x);
  
  // used for unit tests
  static inline bool isEqual(double a, double b) { return (std::abs(a - b) < 1e-5); }
  static unsigned long seed_;
};

} /* namespace maracluster */

#endif /* MARACLUSTER_PVALUECALCULATOR_H_ */
