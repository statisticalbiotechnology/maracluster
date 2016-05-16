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
#include "PvalueCalculator.h"

// to avoid 0 scoring peaks we have to ensure that: 
//
//     log((1-minProb)/minProb) / log((1-maxProb)/maxProb) < PvalueCalculator::probDiscretizationLevels_
//
// for minProb = 1e-10 and maxProb = 0.4 this ratio is 57
unsigned int PvalueCalculator::probDiscretizationLevels_ = 100;
const double PvalueCalculator::kMinProb = 1e-10;
const double PvalueCalculator::kMaxProb = 0.4;
unsigned int PvalueCalculator::kMinScoringPeaks = 15u;
const bool PvalueCalculator::kVariableScoringPeaks = false;

unsigned long PvalueCalculator::seed_ = 1;

void PvalueCalculator::init(const std::vector<unsigned int>& peakBins, const std::vector<double>& peakProbs) {
  peakBins_ = peakBins;
  peakProbs_ = peakProbs;
}

void PvalueCalculator::initPolyfit(const std::vector<unsigned int>& peakBins, const std::vector<unsigned int>& peakScores, const std::vector<double>& polyfit) {
  peakBins_ = peakBins;
  peakScores_ = peakScores;
  maxScore_ = std::accumulate(peakScores_.begin(), peakScores_.end(), 0u);
  polyfit_ = polyfit;
}

void PvalueCalculator::initFromPeakBins(
    const std::vector<unsigned int>& originalPeakBins, 
    const std::vector<double>& peakDist) {
  peakProbs_.clear();
  peakBins_.clear();
  
  BOOST_FOREACH (const unsigned int mzBin, originalPeakBins) {
    if (mzBin >= peakDist.size()) break;
    double peakProb = peakDist.at(mzBin);
    if (peakProb > kMinProb && peakProb < kMaxProb) {
      peakProbs_.push_back(peakProb);
      peakBins_.push_back(mzBin);
      //std::cerr << mzBin << " " << peakProb << std::endl;
    } else if (peakProb <= kMinProb) {
      //std::cerr << "Warning: peak probability below minimum" << std::endl;
    } else if (peakProb >= kMaxProb) {
      //std::cerr << "Warning: peak probability above maximum" << std::endl;
    }
  }
  
  if (peakDist.size() == 0) {
    std::cerr << "Warning: empty peak distribution!" << std::endl;
  }
}

void PvalueCalculator::binaryMatchPeakBins(const std::vector<unsigned int>& queryPeakBins, std::vector<bool>& d) {
  size_t candIdx = 0;
  
  for (size_t i = 0; i < peakBins_.size(); ++i) {
    if (candIdx >= queryPeakBins.size()) break;
    while (peakBins_[i] > queryPeakBins[candIdx]) {
      if (++candIdx >= queryPeakBins.size()) break;
    }
    if (candIdx >= queryPeakBins.size()) break;
    if (peakBins_[i] == queryPeakBins[candIdx]) {
      d[i] = true;
      ++candIdx;
    } else {
      d[i] = false;
    }
  }
}

/**
Input: 
p - list of the probabilities of all the peptides of this protein
d - list including 1 for peptides matched to MS1-features, 0 otherwise
smooth_pvalue - True if the pvalues should be smoothed 
steps_per_position - discretization step

Output: estimated p value
*/
void PvalueCalculator::computePvalVector() {
  
	//std::cerr << "Computing pvalue vector" << std::endl;
  // calculate the vector x and the sum of log(pi)
  std::vector<double> x;
  double sumLogP = 0.0;
  BOOST_FOREACH(double pi, peakProbs_) {
    if (pi >= 0.5) {
      throw std::runtime_error("Found a probability >= 0.5, this will result in an error in pvalue calculation.");
    }
    x.push_back( log( (1.0 - pi)/pi ) );
    sumLogP += log(pi);
  }
  peakProbs_.clear();
  
  // discretize each xi to get li and compute the sum of all li 
  peakScores_.clear();
  unsigned int sumL = 0u;
  double k = 0.0;
  if (x.size() > 0) {
    k = *std::max_element(x.begin(), x.end()) / probDiscretizationLevels_;
    BOOST_FOREACH(double xi, x) {
      unsigned int li = static_cast<unsigned int>(round(xi/k));
      peakScores_.push_back(li);
      if (li == 0) {
        std::cerr << "Warning: zero-scoring peak" << std::endl;
      }
      //std::cout << li << " ";
      sumL += li;
    }
    //std::cout << sumL << " " << l.size() << std::endl;
  }
  maxScore_ = sumL;
  
  // dynamic programming, the easy way
  std::vector<double> f(sumL + 1);
  f[0] = 1.0;
  double c = 0.0;
  unsigned int partSumL = 0u;
  std::vector<unsigned int> ls(peakScores_);
  std::sort(ls.begin(), ls.end());
  for (unsigned int j = 0; j < ls.size(); ++j) {
    if (j % 30 == 29) {
      double fm = *std::max_element(f.begin(), f.begin() + partSumL + 1);
      c += log(fm);
      for (unsigned int i = 0; i <= partSumL; ++i) {
        f[i] /= fm;
      }
    }
    if (ls[j] > 0) {
      partSumL += ls[j];
      for (unsigned int i = partSumL; i >= ls[j]; --i) {
        f[i] += f[i-ls[j]];
      }
    }
  }
  
  // calculate the final p-value vector
  sumProb_.resize(sumL + 1);
  sumProb_[0] = exp(log(f[0]) + sumLogP + c);
  for (unsigned int i = 1; i < sumL + 1; ++i) {
    if (f[i] == 0) {
      sumProb_[i] = sumProb_[i-1];
    } else {
      sumProb_[i] = sumProb_[i-1] + exp (i * k + log(f[i]) + sumLogP + c);
    }
  }
  //std::cerr << "Computed pvalue vector" << std::endl;
}

/**
Input: 
sum_prob - list of the discretized pvalues
d - list with 1s for peak matches, 0s otherwise
l - list of discretized peak logit probabilities
smooth_pvalue - True if the pvalues should be smoothed 

Output: estimated p value
*/
double PvalueCalculator::computePval(const std::vector<unsigned int>& queryPeakBins, bool smoothing) {
  std::vector<bool> d(peakBins_.size());
  binaryMatchPeakBins(queryPeakBins, d);
  
  // express P(D|R = 0) in terms of li (D = obs config)
  double sumThresh = 0.0;
  double sumMax = 0.0;
  unsigned int numMatches = 0;
  double sumScores = 0.0;
  for (unsigned int i = 0; i < d.size(); ++i) {
    if (!d[i]) {
      sumThresh += peakScores_[i];
    } else {
      sumScores += peakScores_[i];
      ++numMatches;
    }
    sumMax += peakScores_[i];
  }
  //std::cout << sumScores/(std::max)(1u,numMatches) << '\t';
  
  //std::cout << numMatches << '\t' << queryPeakBins.size() << "\t" << peakBins_.size() << "\t";
  
  if (smoothing) {
    double lastScoreContrib;
    if (sumThresh > 0) {
      lastScoreContrib = sumProb_[sumThresh] - sumProb_[sumThresh-1];
    } else {
      lastScoreContrib = sumProb_[sumThresh];
    }
    return (std::min)(1.0,sumProb_[sumThresh] - lcg_rand_unif()*lastScoreContrib);
  } else {
    return (std::min)(1.0,sumProb_[sumThresh]);
  }
}

// Based on: http://vilipetek.com/2013/10/07/polynomial-fitting-in-c-using-boost/ (25-07-2014)
void PvalueCalculator::computePvalVectorPolyfit() {
  using namespace boost::numeric::ublas;
  
  computePvalVector();
  
	unsigned int maxScore = sumProb_.size();
	matrix<double> oXMatrix( maxScore, kPolyfitDegree + 1 );
	matrix<double> oYMatrix( maxScore, 1 );
	
	// copy score vector
	for ( size_t i = 0; i < maxScore; i++) {
		oYMatrix(i, 0) = log10(sumProb_[i]);
		//std::cerr << sumProb_[i] << " " << log10(sumProb_[i]) << std::endl;
	}
	sumProb_.clear();

	// create the X (Vandermonde) matrix
	for ( size_t score = 0; score < maxScore; ++score) {
		double value = 1.0;
		double relScore = static_cast<double>(score)/maxScore;
		for ( int col = 0; col < kPolyfitDegree + 1; ++col) {
			oXMatrix(score, col) = value;
			value *= relScore;
		}
	}

	// transpose X matrix
	matrix<double> oXtMatrix( trans(oXMatrix) );
	// multiply transposed X matrix with X matrix
	matrix<double> oXtXMatrix( prec_prod(oXtMatrix, oXMatrix) );
	// multiply transposed X matrix with Y matrix
	matrix<double> oXtYMatrix( prec_prod(oXtMatrix, oYMatrix) );

	// lu decomposition
	permutation_matrix<int> pert(oXtXMatrix.size1());
	const std::size_t singular = lu_factorize(oXtXMatrix, pert);
	// must be singular
	BOOST_ASSERT( singular == 0 );

	// backsubstitution
	lu_substitute(oXtXMatrix, pert, oXtYMatrix);

	// copy the result to coeff
	polyfit_.resize(kPolyfitDegree + 1);
	std::copy( oXtYMatrix.data().begin(), oXtYMatrix.data().end(), polyfit_.begin() );
	
	// TODO: check the residuals?
}

/**
Input: 
sum_prob - list of the discretized pvalues
d - list including 1 for peptides matched to MS1-features, 0 otherwise
l - list of discretized peak logit probabilities
smooth_pvalue - True if the pvalues should be smoothed 

Output: estimated p value
*/
double PvalueCalculator::computePvalPolyfit(const std::vector<unsigned int>& queryPeakBins) {
  std::vector<bool> d(peakBins_.size());
  binaryMatchPeakBins(queryPeakBins, d);
                            
  // express P(D|R = 0) in terms of li (D = obs config)
  unsigned int score = 0u;
  for (unsigned int i = 0; i < d.size(); ++i) {
    if (!d[i]) {
      score += peakScores_[i];
    }/* else {
      std::cerr << "[" << peakBins_[i] << " " << peakScores_[i] << "], ";
    }*/
  }
  
  double relScore = static_cast<double>(score)/maxScore_;
  
  //std::cerr << score << "/" << maxScore_ << " " << polyval(relScore) << std::endl;
  
  return polyval(relScore);
}

double PvalueCalculator::polyval(double x) {
  if (polyfit_.size() > 0) {
    // Horner's method
    double y = polyfit_.back();
    for (int i = polyfit_.size() - 2; i >= 0; --i) {
      y = polyfit_[i] + y*x;
    }
    return (std::min)(0.0,y);
  } else {
    return 0.0;
  }
}

void PvalueCalculator::copyPolyfit(short* peakBins, short* peakScores, double* polyfit) {
  peakBins_.resize(kMaxScoringPeaks, 0u);
  peakScores_.resize(kMaxScoringPeaks, 0u);
  std::copy(peakBins_.begin(), peakBins_.begin() + kMaxScoringPeaks, peakBins);
  std::copy(peakScores_.begin(), peakScores_.begin() + kMaxScoringPeaks, peakScores);
  std::copy(polyfit_.begin(), polyfit_.begin() + kPolyfitDegree + 1, polyfit);
}

void PvalueCalculator::serialize(std::string& polyfitString, std::string& peakScorePairsString) {
  if (polyfit_.size() == 0 || peakBins_.size() == 0 || peakScores_.size() == 0) {
    std::cerr << "Warning: empty vectors in pvalue calculation serialization. " 
              << polyfit_.size() << " " << peakBins_.size() << " " << peakScores_.size() << std::endl;
  }
  polyfitString = "";
  BOOST_FOREACH (double coeff, polyfit_) {
    polyfitString += boost::lexical_cast<std::string>(coeff) + " ";
  }
  polyfitString.substr(0, polyfitString.size() - 1); // remove last space
  
  peakScorePairsString = "";
  if (peakScores_.size() != peakBins_.size()) 
    std::cerr << "Error: peak location and peak score vectors are not of same length" << std::endl;
  
  for (unsigned int j = 0; j < peakScores_.size(); ++j) {
    peakScorePairsString += boost::lexical_cast<std::string>(peakBins_[j]) + " " 
                          + boost::lexical_cast<std::string>(peakScores_[j]) + " ";
  }
  peakScorePairsString.substr(0, peakScorePairsString.size() - 1); // remove last space
}

void PvalueCalculator::deserialize(std::string& polyfitString, std::string& peakScorePairsString) {
  polyfit_.clear();
  std::istringstream iss(polyfitString);
  double coeff;
  while (iss >> coeff) {
    polyfit_.push_back(coeff);
  }
  
  peakBins_.clear();
  peakScores_.clear();
  maxScore_ = 0u;
  
  std::istringstream iss2(peakScorePairsString);
  unsigned int peakBin, score;
  while (iss2 >> peakBin >> score) {
    peakBins_.push_back(peakBin);
    peakScores_.push_back(score);
    maxScore_ += score;
  }
  
  if (polyfit_.size() == 0 || peakBins_.size() == 0 || peakScores_.size() == 0) {
    std::cerr << "Warning: empty vectors in pvalue calculation deserialization. " 
              << polyfit_.size() << " " << peakBins_.size() << " " << peakScores_.size() << std::endl;
  }
}

// Park–Miller random number generator
// from wikipedia
unsigned long PvalueCalculator::lcg_rand() {
  //uint64_t
  seed_ = (seed_ * 279470273) % 4294967291;
  return seed_;
}

double PvalueCalculator::lcg_rand_unif() {
  return (double)(lcg_rand() % 100000) / 100000;
}


bool PvalueCalculator::pvalUnitTest() {
  PvalueCalculator pvalCalc;
  
  std::vector<double> p;
  std::vector<unsigned int> d1, d2;
  
  p.push_back(0.2);
  p.push_back(0.3);
  p.push_back(0.4);
  p.push_back(0.2);
  
  d1.push_back(1);
  d1.push_back(3);
  d1.push_back(5);
  d1.push_back(7);
  
  d2.push_back(1);
  d2.push_back(3);
  d2.push_back(8);
  d2.push_back(9);
  
  pvalCalc.peakProbs_ = p;
  pvalCalc.peakBins_ = d1;
  
  pvalCalc.computePvalVector();
  
  double pval = pvalCalc.computePval(d2);
  
  if (isEqual(pval, 0.135674)) {
    return true;
  } else {
    std::cerr << pval << std::endl;
    std::cerr << 0.135674 << std::endl;
    return false;
  }
}

bool PvalueCalculator::pvalUniformUnitTest() {
  setSeed(100);
  PvalueCalculator pvalCalc;
  
  unsigned int numSpectra = 200000u, numBins = 1000u, numPeaks = 100u;
  double prob = 1-pow(1-1.0/numBins,numPeaks);

  std::vector<double> p(numPeaks, prob);

  std::vector<unsigned int> d1, d2;
  
  std::map<unsigned int, bool> hasPeak;
  for (unsigned int i = 0; i < numPeaks; ++i) {
    bool peakAdded = false;
    while (!peakAdded) {  
      //unsigned int bin = rand() % numBins; // causes biases at 0 (higher density) and 1 (lower density)
      unsigned int bin = lcg_rand() % numBins;
      if (!hasPeak[bin]) {
        d1.push_back(bin);
        peakAdded = true;
      }
    }
  }
  std::sort(d1.begin(), d1.end());
  
  pvalCalc.peakProbs_ = p;
  pvalCalc.peakBins_ = d1;
  
  pvalCalc.computePvalVector();
  
  for (unsigned int j = 0; j < numSpectra; ++j) {
    d2.clear();
    /* produces cap shape
    for (unsigned int i = 0; i < numPeaks; ++i) {  
      unsigned int bin = lcg_rand() % numBins;
      d2.push_back(bin);
    }
    */
    unsigned int hits = 0;
    for (unsigned int i = 0; i < numPeaks; ++i) {  
      if (lcg_rand_unif() < prob) {
        d2.push_back(d1[i]);
        ++hits;
      }
    }
    std::cout << hits << '\t';
    std::sort(d2.begin(), d2.end());
    std::cout << pvalCalc.computePval(d2, true) << std::endl;
  }
  
  return true;
}

bool PvalueCalculator::binaryPeakMatchUnitTest() {
  setSeed(100);
  PvalueCalculator pvalCalc;
  
  unsigned int numBins = 1000u, numPeaks = 100u;
  double prob = 1-pow(1-1.0/numBins,numPeaks);

  std::vector<double> p(numPeaks, prob);

  std::vector<unsigned int> d1, d2;
  
  std::map<unsigned int, bool> hasPeak;
  for (unsigned int i = 0; i < numPeaks; ++i) {
    bool peakAdded = false;
    while (!peakAdded) {  
      //unsigned int bin = rand() % numBins; // causes biases at 0 (higher density) and 1 (lower density)
      unsigned int bin = lcg_rand() % numBins;
      if (!hasPeak[bin]) {
        d1.push_back(bin);
        peakAdded = true;
      }
    }
  }
  std::sort(d1.begin(), d1.end());
  
  for (unsigned int i = 0; i < numPeaks; ++i) {
    bool peakAdded = false;
    while (!peakAdded) {  
      //unsigned int bin = rand() % numBins; // causes biases at 0 (higher density) and 1 (lower density)
      unsigned int bin = lcg_rand() % numBins;
      if (!hasPeak[bin]) {
        d2.push_back(bin);
        peakAdded = true;
      }
    }
  }
  std::sort(d2.begin(), d2.end());
  
  pvalCalc.peakProbs_ = p;
  pvalCalc.peakBins_ = d1;
  
  /*
  for (unsigned int j = 0; j < 1e7; ++j) {
    std::vector<bool> b(d1.size());
    pvalCalc.binaryMatchPeakBins2(d2, b);
  }
  */
  
  unsigned int acc = 0;
  
  std::vector<bool> b(d1.size());
  pvalCalc.binaryMatchPeakBins(d2, b);
    
  BOOST_FOREACH(bool a, b) {
    if (a) ++acc;
  }
  
  if (acc == 12) {
    return true;
  } else {
    std::cout << "Matched peaks was " << acc << ", should be 12." << std::endl;
    return false;
  }
}

bool PvalueCalculator::pvalPolyfitUnitTest() {
  setSeed(10);
  PvalueCalculator pvalCalc;
  
  unsigned int numBins = 1000u, numPeaks = 100u;
  double prob = 1-pow(1-1.0/numBins,numPeaks);

  std::vector<double> p(numPeaks, prob);

  std::vector<unsigned int> d1, d2;
  
  std::map<unsigned int, bool> hasPeak;
  for (unsigned int i = 0; i < numPeaks; ++i) {
    bool peakAdded = false;
    while (!peakAdded) {  
      //unsigned int bin = rand() % numBins; // causes biases at 0 (higher density) and 1 (lower density)
      unsigned int bin = lcg_rand() % numBins;
      if (!hasPeak[bin]) {
        d1.push_back(bin);
        pvalCalc.maxScore_ += 100;
        pvalCalc.peakScores_.push_back(100);
        peakAdded = true;
      }
    }
  }
  std::sort(d1.begin(), d1.end());
  
  for (unsigned int i = 0; i < numPeaks; ++i) {
    bool peakAdded = false;
    while (!peakAdded) {  
      //unsigned int bin = rand() % numBins; // causes biases at 0 (higher density) and 1 (lower density)
      unsigned int bin = lcg_rand() % numBins;
      if (!hasPeak[bin]) {
        d2.push_back(bin);
        peakAdded = true;
      }
    }
  }
  std::sort(d2.begin(), d2.end());
  
  pvalCalc.peakProbs_ = p;
  pvalCalc.peakBins_ = d1;
  pvalCalc.polyfit_.push_back(-40.236851088905013);
  pvalCalc.polyfit_.push_back(90.270774017803348);
  pvalCalc.polyfit_.push_back(-117.64110105643428);
  pvalCalc.polyfit_.push_back(154.88646213751957);
  pvalCalc.polyfit_.push_back(-118.82877111418108);
  pvalCalc.polyfit_.push_back(31.617289887393586);
  
  double logPval = pvalCalc.computePvalPolyfit(d2);
  
  if (isEqual(logPval, -0.817751)) {
    return true;
  } else {
    std::cout << "log(pval) was " << logPval << ", should be -0.817751." << std::endl;
    return false;
  }
}
