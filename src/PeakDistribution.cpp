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
 
/**
This file provides functions for calculating pvalues using dynamic programming
Adapted from https://code.google.com/p/in-silico-mass-fingerprinting/
*/
#include "PeakDistribution.h"

namespace maracluster {

unsigned long PeakDistribution::seed = 1;

void PeakDistribution::init(unsigned int numBins) {
  peakDist_.clear();
  peakDist_.resize(numBins);
}

void PeakDistribution::insert(unsigned int bin, double value) {
  peakDist_.at(bin) = value;
}

void PeakDistribution::print() {
  BOOST_FOREACH(double value, peakDist_) {
    std::cerr << value << " ";
  }
  std::cerr << std::endl;
}

void PeakDistribution::rescale(double scalingFactor) {
  BOOST_FOREACH(double & value, peakDist_) {
    value *= scalingFactor;
  }
}

void PeakDistribution::setUniform(double prob) {
  BOOST_FOREACH(double & value, peakDist_) {
    value = prob;
  }
}

// Calculate Bhattacharyya distance
double PeakDistribution::getDistance(PeakDistribution& otherDist) {
  double bcCoeff = 0.0;
  unsigned int idx = 0;
  unsigned int maxIdx = peakDist_.size();
  BOOST_FOREACH (const double value, otherDist.getDistribution()) {
    if (idx >= maxIdx) break;
    bcCoeff += sqrt(peakDist_[idx] * value);
    //std::cerr << peakDist_[idx] << " " << value << std::endl;
    ++idx;
  }
  //std::cerr << std::endl;
  if (bcCoeff >= 0.0 && bcCoeff <= 1.0 + 1e-5) {
    return -1.0*log(bcCoeff);
  } else {
    std::cerr << "Warning (PeakDistribution::getDistance): Bhattacharyya coefficient out of range, returning maximum distance." << std::endl;
    return 1.0;
  }
}

// Park–Miller random number generator
// from wikipedia
unsigned long PeakDistribution::lcg_rand() {
  //uint64_t
  seed = (seed * 279470273) % 4294967291;
  return seed;
}

double PeakDistribution::lcg_rand_unif() {
  return (double)(lcg_rand() % 100000) / 100000;
}

void PeakDistribution::generateRandSpec(std::vector<unsigned int>& peakBins, const unsigned int numQueryPeaks, bool withDuplicates) {
  peakBins.clear();
  std::vector<double> cumPeakDist(peakDist_.size() + 1);
  cumPeakDist[0] = 0.0;
  
  std::partial_sum(peakDist_.begin(), peakDist_.end(), cumPeakDist.begin() + 1);
  
  std::map<unsigned int, bool> hasPeak;
  for (unsigned int i = 0; i < numQueryPeaks; ++i) {
    bool peakAdded = false;
    while (!peakAdded) {
      double uniRand = lcg_rand_unif();
      for (unsigned int j = 0; j < cumPeakDist.size() - 1; ++j) {
        if (uniRand < cumPeakDist[j+1] && uniRand > cumPeakDist[j]) {
          if (!hasPeak[j] || withDuplicates) {
            peakBins.push_back(j);
            hasPeak[j] = true;
            peakAdded = true;
          }
          break;
        }
      }
    }
  }
  std::sort(peakBins.begin(), peakBins.end());
}

void PeakDistribution::generateRandSpecBernoulli(std::vector<unsigned int>& peakBins, const unsigned int numQueryPeaks, std::vector<unsigned int>& targetPeakBins) {
  peakBins.clear();
  
  unsigned int numMatches = 0;
  BOOST_FOREACH (const unsigned int mzBin, targetPeakBins) {
    if (mzBin >= peakDist_.size()) break;
    //double peakProb = calcPeakProb(peakDist_.at(mzBin), 1, numQueryPeaks);
    double peakProb = peakDist_.at(mzBin);
    if (lcg_rand_unif() < peakProb) {
      peakBins.push_back(mzBin);
      ++numMatches;
    }
  }
  //std::cout << numMatches << '\t';
}

/**
Reads in relative m/z peak probabilities. 
NB: can only parse files with constant step size starting from 0.0 for now
*/
void PeakDistribution::readPeakProbabilities(std::string peakProbFN) {
  std::ifstream peakProbStream(peakProbFN.c_str());
  peakDist_.clear();
  if (peakProbStream.is_open()) {
    bool foundStepSize = false;
    std::string line;
    while (getline(peakProbStream, line)) {
      double tmp, prob;
      std::istringstream iss(line);
      if (!foundStepSize) {
        iss >> stepSize_ >> prob;
        if (stepSize_ > 0.0) {
          foundStepSize = true;
        }
      } else {
        iss >> tmp >> prob;
      }
      peakDist_.push_back(prob);
    }
  } else {
    throw std::runtime_error("Could not open peak probability file");
  }
  peakProbStream.close();
}

void PeakDistribution::getPeakProbabilities(const std::vector<unsigned int>& peakBins, const double precMz, 
                                              std::vector<double>& specPeakProbs, const unsigned int numQueryPeaks) {
  specPeakProbs.clear();
  BOOST_FOREACH (const unsigned int mzBin, peakBins) {
    unsigned int fracBin = getFracBin(mzBin, precMz);
    double peakProb = calcPeakProb(peakDist_[fracBin], precMz, numQueryPeaks);
    if (peakProb > 1e-10 && peakProb < 0.49) {
      specPeakProbs.push_back(peakProb);
    }
  }
}

// clusterp results/140520_Pvalue_clustering/single_peak_probabilities_z3.tsv
bool PeakDistribution::peakReadingUnitTest(std::string peakProbFN) {
  
  PeakDistribution peakDist;
  peakDist.readPeakProbabilities(peakProbFN);
  
  return (peakDist.stepSize_ == 0.0005 && isEqual(peakDist.peakDist_[1000], 0.000118004) );
}

// clusterp results/140520_Pvalue_clustering/single_peak_probabilities_z3.tsv
bool PeakDistribution::peakProbUnitTest(std::string peakProbFN) {
  std::vector<double> specPeakProbs;  
  PeakDistribution peakDist;
  
  peakDist.readPeakProbabilities(peakProbFN);
  
  std::vector<unsigned int> peakBins;
  
  peakBins.push_back(200);
  peakBins.push_back(400);
  peakBins.push_back(600);
  peakDist.getPeakProbabilities(peakBins, 500.0, specPeakProbs, 100u);
  
  return (specPeakProbs.size() == 3 && isEqual(specPeakProbs[0], 0.0156541) && isEqual(specPeakProbs[1], 0.0862822) && isEqual(specPeakProbs[2], 0.13269) );
}

} /* namespace maracluster */
