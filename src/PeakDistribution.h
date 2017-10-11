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
 
#ifndef MARACLUSTER_PEAKDISTRIBUTION_H_
#define MARACLUSTER_PEAKDISTRIBUTION_H_

#include <vector>
#include <map>
#include <string>

#include <cstdlib>
#include <stdexcept>

#include <cmath>
#include <algorithm>
#include <numeric>

#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/foreach.hpp>

namespace maracluster {

class PeakDistribution {
 public:    
  void init(unsigned int numBins);
  void insert(unsigned int bin, double value);
  void rescale(double scalingFactor);
  void setUniform(double prob);
  
  const std::vector<double>& getDistribution() const { return peakDist_; }
  double getDistance(PeakDistribution& otherDist);
  
  void print();
  
  void readPeakProbabilities(std::string peakProbFN);
  void getPeakProbabilities(const std::vector<unsigned int>& peakBins, 
    const double precMz, std::vector<double>& specPeakProbs,
    const unsigned int numQueryPeaks);
  
  void getPeakProbabilities(std::vector<unsigned int>& peakBins, 
    std::vector<double>& specPeakProbs, const unsigned int numQueryPeaks);
  
  void generateRandSpec(std::vector<unsigned int>& peakBins, const unsigned int numPeaks = 100, bool withDuplicates = false);
  void generateRandSpecBernoulli(std::vector<unsigned int>& peakBins, const unsigned int numQueryPeaks, 
                                 std::vector<unsigned int>& targetPeakBins);
  
  inline static void setSeed(unsigned long s) { seed = s; }
  static unsigned long lcg_rand();
  static double lcg_rand_unif();
  
  // clusterp results/140520_Pvalue_clustering/single_peak_probabilities_z3.tsv
  static bool peakReadingUnitTest(std::string peakProbFN);
  static bool peakProbUnitTest(std::string peakProbFN);
 private:
  std::vector<double> peakDist_;
  double stepSize_;
  
  inline unsigned int getFracBin(unsigned int mzBin, double precMz) {
    return static_cast<unsigned int>( mzBin / (precMz * stepSize_) );
  }
  
  static inline double calcPeakProb(double singlePeakProb, double precMz = 1, unsigned int nPeaks = 100) {
    return (1 - std::pow(1 - 1/precMz * singlePeakProb, static_cast<int>(nPeaks)));
  }
  
  // used for unit tests
  static inline bool isEqual(double a, double b) { return (std::abs(a - b) < 1e-5); }
  
  static unsigned long seed;
};

} /* namespace maracluster */

#endif /* MARACLUSTER_PEAKDISTRIBUTION_H_ */
