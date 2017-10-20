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
 
#include "BatchPvalueVectors.cuh"

namespace maracluster {

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__device__
void polyval(double *polyfit, double *x) {
  // Horner's method
  double y = polyfit[POLYFIT_SIZE - 1];
  for (int i = POLYFIT_SIZE - 2; i >= 0; --i) {
    y = polyfit[i] + y*(*x);
  }
  if (y > 0.0) {
    *x = 0.0;
  } else {
    *x = y;
  }
}

__device__
void binaryMatchPeakBins(short *targetPeakBins, short *targetPeakScores, short *queryPeakBins, int *score) {
  size_t qIdx = 0;
  for (size_t i = 0; i < SCORING_PEAKS; ++i) {
    if (qIdx >= SCORING_PEAKS) break;
    while (targetPeakBins[i] > queryPeakBins[qIdx]) {
      if (++qIdx >= SCORING_PEAKS) break;
    }
    if (qIdx >= SCORING_PEAKS) break;
    if (targetPeakBins[i] == queryPeakBins[qIdx]) {
      ++qIdx;
    } else {
      *score += targetPeakScores[i];
    }
  }
}

__global__
void calculatePvals(short *peakBins, short *peakScores, double *polyfits, int *maxScores, int n, short *queryPeakBins, int m, double *pvals) {
  size_t targetIdx = blockIdx.x * BLOCK_SIZE + threadIdx.x;
  size_t queryIdx = blockIdx.y * BLOCK_SIZE + threadIdx.y;
  
  if (targetIdx < n && queryIdx < m) {
    size_t flattenedTargetIdx = targetIdx * SCORING_PEAKS;
    size_t flattenedQueryIdx = queryIdx * SCORING_PEAKS;
    
    int score = 0;
    binaryMatchPeakBins(&peakBins[flattenedTargetIdx], &peakScores[flattenedTargetIdx], &queryPeakBins[flattenedQueryIdx], &score);
    
    double relScore = static_cast<double>(score) / maxScores[targetIdx];
    polyval(&polyfits[targetIdx * POLYFIT_SIZE], &relScore);
    size_t pvalIdx = targetIdx * m + queryIdx;
    pvals[pvalIdx] = relScore;
  }
}

__host__
void runKernel(short *peakBins, short *peakScores, double *polyfits, int *maxScores, size_t N, short *queryPeakBins, size_t M, double *pvals) {
  short *peakBinsDevice, *peakScoresDevice, *queryPeakBinsDevice;
  double *polyfitsDevice, *pvalsDevice;
  int *maxScoresDevice;
  
  gpuErrchk( cudaMalloc(&peakBinsDevice, N * SCORING_PEAKS * sizeof(short)) );
  gpuErrchk( cudaMalloc(&peakScoresDevice, N * SCORING_PEAKS * sizeof(short)) );
  gpuErrchk( cudaMalloc(&polyfitsDevice, N * POLYFIT_SIZE * sizeof(double)) );
  gpuErrchk( cudaMalloc(&maxScoresDevice, N * sizeof(int)) );
  gpuErrchk( cudaMalloc(&pvalsDevice, N * PVEC_MAX_BATCH_SIZE * sizeof(double)) );
  gpuErrchk( cudaMalloc(&queryPeakBinsDevice, M * SCORING_PEAKS * sizeof(short)) );
  
  gpuErrchk( cudaMemcpy(peakBinsDevice, peakBins, N * SCORING_PEAKS * sizeof(short), cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(peakScoresDevice, peakScores, N * SCORING_PEAKS * sizeof(short), cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(polyfitsDevice, polyfits, N * POLYFIT_SIZE * sizeof(double), cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(maxScoresDevice, maxScores, N * sizeof(int), cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(queryPeakBinsDevice, queryPeakBins, M * SCORING_PEAKS * sizeof(short), cudaMemcpyHostToDevice) );
  
  dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
  dim3 dimGrid;
  dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
  dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;
  
  calculatePvals<<<dimGrid, dimBlock>>>(peakBinsDevice, peakScoresDevice, polyfitsDevice, maxScoresDevice, N, queryPeakBinsDevice, M, pvalsDevice);
  
  //gpuErrchk( cudaDeviceSynchronize() );
  
  gpuErrchk( cudaMemcpy(pvals, pvalsDevice, N * PVEC_MAX_BATCH_SIZE * sizeof(double), cudaMemcpyDeviceToHost) );
  
  /*
  double sumPvals = 0;
  for (size_t i = 0; i < N; i++) {
    sumPvals += pvals[i*M];
  }
  std::cout << "Sum: " << sumPvals << std::endl;
  */
  
  // Free memory
  gpuErrchk( cudaFree(peakBinsDevice) );
  gpuErrchk( cudaFree(peakScoresDevice) );
  gpuErrchk( cudaFree(polyfitsDevice) );
  gpuErrchk( cudaFree(maxScoresDevice) );
  gpuErrchk( cudaFree(pvalsDevice) );
  gpuErrchk( cudaFree(queryPeakBinsDevice) );
}

} /* namespace maracluster */
