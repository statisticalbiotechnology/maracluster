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
 
#include "ClusterMerge.h"

namespace maracluster {

double ClusterMerge::mzWeight = 1.0;
double ClusterMerge::intWeight = 0.1;

double ClusterMerge::peakDistance(MZIntensityPair mzi1, MZIntensityPair mzi2) {
  return std::abs(mzi1.mz - mzi2.mz)*mzWeight + std::abs(mzi1.intensity - mzi2.intensity)*intWeight;
}

MZIntensityPair ClusterMerge::mergePeaks(MZIntensityPair mzi1, MZIntensityPair mzi2) {
  return MZIntensityPair( 
    (mzi1.multiplicity*mzi1.mz + mzi2.multiplicity*mzi2.mz)/(mzi1.multiplicity+mzi2.multiplicity), 
    mzi1.intensity + mzi2.intensity, 
    mzi1.multiplicity+mzi2.multiplicity 
  );
}
  
MZIntensityPair ClusterMerge::findNN(MZIntensityPair mzi1, const std::vector<MZIntensityPair>& mziPairsFrom) {
  MZIntensityPair nnMzi(0.0,0.0);
  double nnDist = 1000000.0;
  BOOST_FOREACH (const MZIntensityPair& mzi2, mziPairsFrom) {
    double dist = peakDistance(mzi1, mzi2);
    if (dist < nnDist) {
      nnMzi = mzi2;
      nnDist = dist;
    }
  }
  return nnMzi;
}
  
void ClusterMerge::merge(std::vector<MZIntensityPair>& mziPairsIn, std::vector<MZIntensityPair>& mziPairsFrom, double weight) {
  std::vector<MZIntensityPair> mziPairsMerged(mziPairsIn);
  double sumAbsIntDiff = 0.0;
  double sumAbsMzDiff = 0.0;
  BOOST_FOREACH (MZIntensityPair& mzi2, mziPairsFrom) {
    MZIntensityPair nnMzi = findNN(mzi2, mziPairsIn);
    if (std::abs(nnMzi.mz - mzi2.mz) <= 0.3 && std::abs(nnMzi.intensity/nnMzi.multiplicity - mzi2.intensity/mzi2.multiplicity) <= 5.0) {
      std::cerr << nnMzi.mz - mzi2.mz << " ";
      std::cerr << nnMzi.intensity/nnMzi.multiplicity - mzi2.intensity/mzi2.multiplicity << std::endl;
      sumAbsIntDiff += std::abs(nnMzi.intensity - mzi2.intensity);
      sumAbsMzDiff += std::abs(nnMzi.mz - mzi2.mz);
      MZIntensityPair mergedMzi = mergePeaks(mzi2, nnMzi);
      mziPairsMerged.erase(std::remove(mziPairsMerged.begin(), mziPairsMerged.end(), nnMzi), mziPairsMerged.end());
      //print "\t".join(map(str,[mz1, int1, nnMz, nnInt, mergedMz, mergedInt]))
      mziPairsMerged.push_back(mergedMzi);
    } else {
      mziPairsMerged.push_back(mzi2);
    }
  }
  //std::cerr << mziPairsMerged.size() << std::endl;
  //std::cerr << "Avg intensity diff = " << sumAbsIntDiff / mziPairsMerged.size() << "; avg mz diff = " << sumAbsMzDiff / mziPairsMerged.size() << " Number of merged peaks = " << mziPairsMerged.size() << std::endl;
  std::sort( mziPairsMerged.begin(), mziPairsMerged.end(), SpectrumHandler::greaterMultiplicity );
  mziPairsIn.clear();
  std::copy( mziPairsMerged.begin(), mziPairsMerged.begin() + (std::min)(200, static_cast<int>(mziPairsMerged.size())), std::back_inserter(mziPairsIn) );
  //std::cerr << mziPairsIn.size() << " " << *mziPairsIn.begin() << std::endl;
}
  

} /* namespace maracluster */
