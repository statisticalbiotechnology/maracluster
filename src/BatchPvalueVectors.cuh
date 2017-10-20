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

#ifndef MARACLUSTER_BATCHPVALUEVECTORS_CUH_
#define MARACLUSTER_BATCHPVALUEVECTORS_CUH_

/* doubling this batch size from 16384 to 32768 runs out of memory */
#define PVEC_MAX_BATCH_SIZE 1024
#define BLOCK_SIZE 16
#define SCORING_PEAKS 40
#define POLYFIT_SIZE 6

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdio>

namespace maracluster {

// Kernel function to add the elements of two arrays
void runKernel(short *peakBins, short *peakScores, double *polyfits, int *maxScores, size_t N, short *queryPeakBins, size_t M, double* pvals);

} /* namespace maracluster */

#endif /* MARACLUSTER_BATCHPVALUEVECTORS_CUH_ */
