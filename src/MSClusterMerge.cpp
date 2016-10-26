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
 
#include "MSClusterMerge.h"

const float MSClusterMerge::MASS_TO_INT_RATIO = 10000.0f;
float MSClusterMerge::fragmentTolerance_ = 0.34f; /* has to be > 0.1 */
float MSClusterMerge::isoTolerance_ = 
    0.1f + (MSClusterMerge::fragmentTolerance_ - 0.1f) * 0.5f;
PeakWeightTable MSClusterMerge::peakWeightTable_;

void MSClusterMerge::init() {
  peakWeightTable_.initWeights(64, 0.15f);
}

/* Adapted from the MS-Cluster code base:
   software/tools/ms-cluster/src/PepNovo/PeakList.cpp joinAdjacentPeaks() */
void MSClusterMerge::binMZIntensityPairs(std::vector<MZIntensityPair>& spectrum,
    std::vector<BinnedMZIntensityPair>& binnedSpectrum) {
  binnedSpectrum.clear();
  if (spectrum.size() == 0) return;
  
  binnedSpectrum.reserve(spectrum.size());
  
  int totalPeaks = static_cast<int>(spectrum.size());
	const int maxIntProximity = convertMassToInt(fragmentTolerance_);
	binnedSpectrum.push_back(BinnedMZIntensityPair(
	    convertMassToInt(spectrum[0].mz), spectrum[0].mz, spectrum[0].intensity));
	int prev = 0;
  for (int i = 1; i < totalPeaks; i++) {
    int curBinIdx = convertMassToInt(spectrum[i].mz);
		if (curBinIdx - static_cast<int>(binnedSpectrum[prev].binIdx) < maxIntProximity) {
			// join peaks with proportion to their intensities
			const double intensitySum = 
			    (binnedSpectrum[prev].intensity + spectrum[i].intensity);
			const double ratio = binnedSpectrum[prev].intensity/intensitySum;
			const double newMass = ratio * binnedSpectrum[prev].mz + 
                             (1.0-ratio) * spectrum[i].mz;
			
			binnedSpectrum[prev].mz = newMass;
			binnedSpectrum[prev].binIdx = convertMassToInt(newMass);
			binnedSpectrum[prev].intensity = intensitySum;
		} else {
			binnedSpectrum.push_back(BinnedMZIntensityPair(
			    curBinIdx, spectrum[i].mz, spectrum[i].intensity));
			++prev;
		}
	}
  
  binnedSpectrum.resize(prev+1);
}

void MSClusterMerge::unbinMZIntensityPairs(
    std::vector<BinnedMZIntensityPair>& binnedSpectrum,
    std::vector<MZIntensityPair>& spectrum) {
  spectrum.clear();
  
  BOOST_FOREACH (const MZIntensityPair& mziPair, binnedSpectrum) {
    spectrum.push_back(MZIntensityPair(mziPair.mz, mziPair.intensity));
  }
  
}

void MSClusterMerge::mergeMccs(std::vector<MassChargeCandidate>& allMccs,
    std::vector<MassChargeCandidate>& consensusMccs) {  
  std::sort(allMccs.begin(), allMccs.end(), MassChargeCandidate::lessChargeMass);
  
  std::vector<unsigned int> mccCounts;
  std::vector<MassChargeCandidate> candidateMccs;
  candidateMccs.push_back(allMccs.front());
  candidateMccs.back().precMz = SpectrumHandler::calcPrecMz(
		    candidateMccs.back().mass,
        candidateMccs.back().charge);
  mccCounts.push_back(1);
  for (int i = 1; i < static_cast<int>(allMccs.size()); ++i) {
    double massTolerance = 10e-6 * candidateMccs.back().mass; /* 10 ppm */
    if (allMccs[i].charge == candidateMccs.back().charge &&
        allMccs[i].mass - candidateMccs.back().mass < massTolerance) {
			// join peaks with proportion to their multiplicities
			const double ratio = 1.0/mccCounts.back();
			const double newMass = ratio * allMccs[i].mass + 
                             (1.0-ratio) * candidateMccs.back().mass;

			candidateMccs.back().mass = newMass;
			mccCounts.back() += 1;
		} else {
			candidateMccs.push_back(allMccs[i]);
			mccCounts.push_back(1);
		}
		candidateMccs.back().precMz = SpectrumHandler::calcPrecMz(
		    candidateMccs.back().mass,
        candidateMccs.back().charge);
  }
  
  consensusMccs.clear();  
  //consensusMccs.swap(candidateMccs);
  
  // keep only the top 3 most frequent mccs, keeping ties for third
  std::vector<unsigned int> mccCountsCopy = mccCounts;
  std::sort(mccCountsCopy.begin(), mccCountsCopy.end());
  unsigned int threshold = mccCountsCopy[mccCountsCopy.size()-std::min(mccCountsCopy.size(),static_cast<size_t>(3))];
  //unsigned int maxCount = mccCountsCopy[mccCountsCopy.size()-1];
  
  for (size_t i = 0; i < candidateMccs.size(); ++i) {
    if (mccCounts[i] >= threshold/* && mccCounts[i]*3 >= maxCount */) {
      consensusMccs.push_back(candidateMccs[i]);
    }
  }
  
}

/* Adapted from the MS-Cluster code base:
   software/tools/ms-cluster/src/MsCluster/MsClusterDataStorage.cpp */
void MSClusterMerge::merge(
    std::vector< std::vector<BinnedMZIntensityPair> >& cluster,
    std::vector<BinnedMZIntensityPair>& mergedSpectrum) {
	const size_t peakAreaSize = 160;
	//const size_t peakAreaSize = 1000;
	std::vector<BinnedMZIntensityPair> allPeaks;
	std::vector<int> peakCounts;
	
	int totalPeaks = 0;
	for (size_t i = 0; i < cluster.size(); ++i) {
		int numPeaks = 0;
		for (size_t j = 0; j < cluster[i].size(); ++j) {
		  if (cluster[i][j].intensity > 0) {
		    allPeaks.push_back(cluster[i][j]);
		    peakCounts.push_back(1);
		    ++numPeaks;
		  }
	  }
		totalPeaks += numPeaks;
	}
  
  // sort by mass
	sort(allPeaks.begin(), allPeaks.end(), SpectrumHandler::lessMZ);

	// merge peaks
	mergedSpectrum.clear();
	if (allPeaks.size() == 0) return;
	
	mergedSpectrum.push_back(allPeaks.at(0));
	const int maxIntProximity = convertMassToInt(isoTolerance_);
	int prev = 0;
	for (int i = 1; i < totalPeaks; i++) {
		if (static_cast<int>(allPeaks.at(i).binIdx - mergedSpectrum.at(prev).binIdx) < maxIntProximity ) {
			// join peaks with proportion to their intensities
			const double intensitySum = 
			    (mergedSpectrum.at(prev).intensity + allPeaks.at(i).intensity);
			const double ratio = mergedSpectrum.at(prev).intensity/intensitySum;
			const double newMass = ratio * mergedSpectrum.at(prev).mz + 
                             (1.0-ratio) * allPeaks.at(i).mz;
			
			mergedSpectrum.at(prev).mz = newMass;
			mergedSpectrum.at(prev).binIdx = convertMassToInt(newMass);
			mergedSpectrum.at(prev).intensity = intensitySum;
			peakCounts.at(prev) += peakCounts.at(i);
		} else {
			mergedSpectrum.push_back(allPeaks.at(i));
			peakCounts.at(++prev) = peakCounts.at(i);
		}
	}
	totalPeaks = prev + 1;
  
	// modify the intensity according to the peakWeightTable_
	// that is discount the weight of peaks that have only a few copies
	for (int i = 0; i < totalPeaks; i++) {
	  //std::cerr << mergedSpectrum[i].mz << " " <<  mergedSpectrum[i].intensity << " " << peakCounts[i] << " " << peakWeightTable_.getWeight(static_cast<int>(peakCounts[i]), static_cast<int>(cluster.size())) << std::endl;
		mergedSpectrum.at(i).intensity *=
			  peakWeightTable_.getWeight(static_cast<int>(peakCounts.at(i)),
			                             static_cast<int>(cluster.size()));
  }
	// select a number of peaks according to their intensity
	sort(mergedSpectrum.begin(), mergedSpectrum.end(),
	     SpectrumHandler::greaterIntensity);

  mergedSpectrum.resize((std::min)(peakAreaSize, mergedSpectrum.size()));
  
	// sort according to mass
	sort(mergedSpectrum.begin(), mergedSpectrum.end(), SpectrumHandler::lessMZ);
}

bool MSClusterMerge::mergeUnitTest() {
  std::vector<MZIntensityPair> mziPairs1, mziPairs2, mziPairs3;
  std::vector<BinnedMZIntensityPair> spectrum1, spectrum2, spectrum3;
  
  mziPairs1.push_back(MZIntensityPair(1000,10.0));
  mziPairs1.push_back(MZIntensityPair(2000,15.0));
  mziPairs1.push_back(MZIntensityPair(3000,20.0));
  binMZIntensityPairs(mziPairs1, spectrum1);
  
  mziPairs2.push_back(MZIntensityPair(3000.01,10.0));
  mziPairs2.push_back(MZIntensityPair(4000,15.0));
  mziPairs2.push_back(MZIntensityPair(5000.01,20.0));
  binMZIntensityPairs(mziPairs2, spectrum2);
  
  mziPairs3.push_back(MZIntensityPair(5000,10.0));
  mziPairs3.push_back(MZIntensityPair(6000,15.0));
  mziPairs3.push_back(MZIntensityPair(7000,20.0));
  binMZIntensityPairs(mziPairs3, spectrum3);
  
  std::vector< std::vector<BinnedMZIntensityPair> > cluster;
  std::vector<BinnedMZIntensityPair> mergedSpectrum;
  cluster.push_back(spectrum1);
  cluster.push_back(spectrum2);
  cluster.push_back(spectrum3);
  
  init();
  merge(cluster, mergedSpectrum);
  
  /*
  10000000 6.14125
  20000000 9.21188
  30000034 30
  40000000 9.21188
  50000068 30
  60000000 9.21188
  70000000 12.2825
  */
  
  unsigned int numFound = 0u;
  bool wrongIntensity = false;
  BOOST_FOREACH(BinnedMZIntensityPair& p, mergedSpectrum) {
    if (p.binIdx == 10000000) {
      ++numFound;
      if (!isEqual(p.intensity, 6.14125)) {
        std::cerr << "Intensity for binIdx 10000000 is " << p.intensity << ", should be 6.14125." << std::endl;
        wrongIntensity = true;
      }
    } else if (p.binIdx == 30000034) {
      ++numFound;
      if (!isEqual(p.intensity, 30.0)) {
        std::cerr << "Intensity for binIdx 30000034 is " << p.intensity << ", should be 30.0." << std::endl;
        wrongIntensity = true;
      }
    } else if (p.binIdx == 50000068) {
      ++numFound;
      if (!isEqual(p.intensity, 30.0)) {
        std::cerr << "Intensity for binIdx 50000068 is " << p.intensity << ", should be 30.0." << std::endl;
        wrongIntensity = true;
      }
    }
	  //std::cout << p.binIdx << " " << p.intensity << std::endl;
	}
	
	if (wrongIntensity) return false;
	else if (numFound < 3) {
	  std::cerr << "One of the peaks is missing in the consensus spectrum." << std::endl;
	  return false;
	}
	
  return true;
}
