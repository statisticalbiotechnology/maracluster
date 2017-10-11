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
 
#include "InterpolationMerge.h"

namespace maracluster {

//double InterpolationMerge::maxScaling = 0.001;
//double InterpolationMerge::maxShift = 0.5;
double InterpolationMerge::maxScaling = 0.0;
double InterpolationMerge::maxShift = 0.0;

void InterpolationMerge::merge(std::vector<MZIntensityPair>& mziPairsIn, std::vector<MZIntensityPair>& mziPairsFrom, 
                                  double weight, Method method, double shift, double scaling) {
  unsigned int fromCount = 0;
  BOOST_FOREACH (MZIntensityPair& mziPair, mziPairsIn) {
    while (mziPairsFrom.at(fromCount).mz < (mziPair.mz + shift)*scaling) {
      ++fromCount;
      if (fromCount >= mziPairsFrom.size()) break;
    }
    if (fromCount >= mziPairsFrom.size()) break;
    if (fromCount >= 1) {
      switch (method) {
        case MINM:
        {
          double interpolatedIntensity = SpectrumHandler::interpolateIntensity(mziPairsFrom.at(fromCount-1), mziPairsFrom.at(fromCount), (mziPair.mz + shift)*scaling);
          mziPair.intensity = (std::min)(interpolatedIntensity, mziPair.intensity);
          break;
        }
        case MAXM:
        {
          double interpolatedIntensity = SpectrumHandler::interpolateIntensity(mziPairsFrom.at(fromCount-1), mziPairsFrom.at(fromCount), (mziPair.mz + shift)*scaling);
          mziPair.intensity = (std::max)(interpolatedIntensity, mziPair.intensity);
          break;
        }
        default: // AVGM
        {
          mziPair.intensity += weight * SpectrumHandler::interpolateIntensity(mziPairsFrom.at(fromCount-1), mziPairsFrom.at(fromCount), (mziPair.mz + shift)*scaling);
          break;
        }
      }
    }
  }
}

void InterpolationMerge::mergeViceVersa(std::vector<MZIntensityPair>& mziPairsIn, std::vector<MZIntensityPair>& mziPairsFrom, 
                                            double weight, double shift, double scaling) {
  std::vector<MZIntensityPair> mziPairsInCopy = mziPairsIn;
  merge(mziPairsIn, mziPairsFrom, weight, AVGM, shift, scaling);
  merge(mziPairsFrom, mziPairsInCopy, 1.0/weight, AVGM, -shift, 1.0/scaling);
  SpectrumHandler::scaleIntensities(mziPairsFrom, weight);
  mziPairsIn.reserve( mziPairsIn.size() + mziPairsFrom.size() );
  mziPairsIn.insert( mziPairsIn.end(), mziPairsFrom.begin(), mziPairsFrom.end() );
}

double InterpolationMerge::interpolateDotProduct(const std::vector<MZIntensityPair>& mziPairsIn, const std::vector<MZIntensityPair>& mziPairsFrom, double weight, double shift, double scaling) {
  unsigned int fromCount = 0;
  double dotProduct = 0.0;
  BOOST_FOREACH (MZIntensityPair mziPair, mziPairsIn) {
    while (mziPairsFrom.at(fromCount).mz < (mziPair.mz + shift)*scaling) {
      ++fromCount;
      if (fromCount >= mziPairsFrom.size()) break;
    }
    if (fromCount >= mziPairsFrom.size()) break;
    if (fromCount >= 1) {
      double prod = mziPair.intensity * SpectrumHandler::interpolateIntensity(mziPairsFrom.at(fromCount-1), mziPairsFrom.at(fromCount), (mziPair.mz + shift)*scaling);
      dotProduct += prod * prod;
    }
  }
  return dotProduct;
}

void InterpolationMerge::findBestAffineTransform(std::vector<MZIntensityPair>& mziPairsIn, std::vector<MZIntensityPair>& mziPairsFrom, 
                                  const double weight, double & winningShift, double & winningScaling) {
  winningShift = 0.0;
  winningScaling = 1.0;
  
  if (maxScaling != 0 || maxShift != 0.0) {
    double winningDotProduct = 0.0;
    for (double scaling = 1-maxScaling; scaling <= 1+maxScaling; scaling += maxScaling/50) {
      for (double shift = -maxShift; shift <= maxShift; shift += maxShift/50) {
        double dotProduct = interpolateDotProduct(mziPairsIn, mziPairsFrom, weight, shift, scaling);
        if (dotProduct > winningDotProduct) {
          winningDotProduct = dotProduct;
          winningShift = shift;
          winningScaling = scaling;
        }
      }
    }
  }
}

} /* namespace maracluster */
