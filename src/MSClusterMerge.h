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
 
#ifndef MARACLUSTER_MSCLUSTERMERGE_H_
#define MARACLUSTER_MSCLUSTERMERGE_H_

#include <vector>
#include <cmath>

#include <boost/foreach.hpp>

#include "SpectrumHandler.h"
#include "BinSpectra.h"

namespace maracluster {

/* Adapted from software/tools/ms-cluster/src/MsCluster/PeakWeightTable.h */
class PeakWeightTable {
 public:
	void initWeights(int maxN, float p) {
		assert(p>0.0 && p<1.0);
		maxN_ = maxN;
		kWithValueOne_.resize(maxN_+1);
		weights_.resize(maxN_+1);
		for (int n=1; n<=maxN_; n++) {
			weights_[n].push_back(0.0);
			std::vector<float> cdf;
			computeBinomialCDFs(n, p, cdf);
			int i;
			for (i=0; i<=n; i++) {
				float w=cdf[i];
				if (w > 0.99)
					break;
				weights_[n].push_back(w);
			}
			kWithValueOne_[n]=i;
		}
	}

	float getWeight(int k, int n) {
		while (n>maxN_) {
			n >>= 1;
			k >>= 1;
		}

		if (k>=kWithValueOne_[n])
			return 1.0;
		return (weights_[n][k]);
	}

 private:
	int maxN_;
	std::vector<int> kWithValueOne_;
	/* the weights for having k of n peaks */
	std::vector< std::vector<float> > weights_; 
	
	void computeBinomialCDFs(int n, float p, std::vector<float>& cdf) {
	  cdf.clear();
	  if (n<=0 || p<=0.0 || p>=1.0)
		  return;

	  cdf.resize(n+1,0.0);
	  const double pOverOneMinusP = p/(1.0-p);
	  double binomialCoeff = 1.0;				// n choose k
	  double powerValue    = pow(1.0-p, n);   // p^k * (1-p)^n-k
	  cdf[0] = static_cast<float>(powerValue);

	  for (int k=1; k<=n; k++) {
		  binomialCoeff *= static_cast<double>(n-k+1.0);
		  binomialCoeff /= static_cast<double>(k);
		  powerValue	  *= pOverOneMinusP;
		  cdf[k]= cdf[k-1] + static_cast<float>(binomialCoeff * powerValue);
	  }
  }
};

class MSClusterMerge {
 public:
  static float fragmentTolerance_;
  static float isoTolerance_;
  
  static void init();
  static void merge(std::vector< std::vector<BinnedMZIntensityPair> >& cluster,
                    std::vector<BinnedMZIntensityPair>& mergedSpectrum);
  
  static void mergeMccs(std::vector<MassChargeCandidate>& allMccs,
      std::vector<MassChargeCandidate>& consensusMccs);
  static bool mergeUnitTest();
  
  static void binMZIntensityPairs(std::vector<MZIntensityPair>& spectrum,
      std::vector<BinnedMZIntensityPair>& binnedSpectrum);
  static void unbinMZIntensityPairs(
      std::vector<BinnedMZIntensityPair>& binnedSpectrum,
      std::vector<MZIntensityPair>& spectrum);
 private:
  static PeakWeightTable peakWeightTable_; // for creating consensuses
  
  static const float MASS_TO_INT_RATIO;

  inline static int convertMassToInt(float mass) {
	  return (static_cast<int>(MASS_TO_INT_RATIO * mass));
  }
  
  inline static int convertMassToInt(double mass) {
	  return convertMassToInt(static_cast<float>(mass));
  }

  inline static float convertIntToMass(int m) {
	  return (static_cast<float>(m)/MASS_TO_INT_RATIO);
  }
  
  inline static bool lessMZ(const std::pair<BinnedMZIntensityPair, size_t>& a, 
                            const std::pair<BinnedMZIntensityPair, size_t>& b) { 
    return SpectrumHandler::lessMZ(a.first, b.first); 
  }
  
  // used for unit tests
  static inline bool isEqual(double a, double b) { return (std::abs(a - b) < 1e-5); }
};

} /* namespace maracluster */

#endif /* MARACLUSTER_MSCLUSTERMERGE_H_ */
