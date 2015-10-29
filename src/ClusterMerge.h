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
 
#ifndef CLUSTER_MERGE_H
#define CLUSTER_MERGE_H

#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <cmath>

#include <boost/foreach.hpp>

#include "SpectrumHandler.h"
#include "MZIntensityPair.h"

class ClusterMerge {
  public:
    static void merge(std::vector<MZIntensityPair>& mziPairsIn, std::vector<MZIntensityPair>& mziPairsFrom, double weight = 1.0);
  protected:
    static double mzWeight, intWeight;
    static double peakDistance(MZIntensityPair mzi1, MZIntensityPair mzi2);
    static MZIntensityPair mergePeaks(MZIntensityPair mzi1, MZIntensityPair mzi2);
    static MZIntensityPair findNN(MZIntensityPair mzi1, const std::vector<MZIntensityPair>& mziPairsFrom);
};

#endif // CLUSTER_MERGE_H
