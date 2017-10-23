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

float totaltime = 0.0f;
cudaStream_t streams[NUM_STREAMS];
short *peakBinsDevice[NUM_STREAMS];
short *peakScoresDevice[NUM_STREAMS];
short *queryPeakBinsDevice[NUM_STREAMS];
double *polyfitsDevice[NUM_STREAMS];
double *pvalsDevice[NUM_STREAMS];
int *maxScoresDevice[NUM_STREAMS];
  
void initStreams(double **pvalsHost) {
  for (int i = 0; i < NUM_STREAMS; ++i) { 
    gpuErrchk( cudaStreamCreate(&streams[i]) );
    gpuErrchk( cudaMalloc(&peakBinsDevice[i], PVEC_MAX_BATCH_SIZE * SCORING_PEAKS * sizeof(short)) );
    gpuErrchk( cudaMalloc(&peakScoresDevice[i], PVEC_MAX_BATCH_SIZE * SCORING_PEAKS * sizeof(short)) );
    gpuErrchk( cudaMalloc(&polyfitsDevice[i], PVEC_MAX_BATCH_SIZE * POLYFIT_SIZE * sizeof(double)) );
    gpuErrchk( cudaMalloc(&maxScoresDevice[i], PVEC_MAX_BATCH_SIZE * sizeof(int)) );
    gpuErrchk( cudaMalloc(&pvalsDevice[i], PVEC_MAX_BATCH_SIZE * PVEC_MAX_BATCH_SIZE * sizeof(double)) );
    gpuErrchk( cudaMallocHost(&pvalsHost[i], PVEC_MAX_BATCH_SIZE * PVEC_MAX_BATCH_SIZE * sizeof(double)) );
    gpuErrchk( cudaMalloc(&queryPeakBinsDevice[i], PVEC_MAX_BATCH_SIZE * SCORING_PEAKS * sizeof(short)) );
  }
}

void destroyStreams(double **pvalsHost) {
  for (int i = 0; i < NUM_STREAMS; ++i) { 
    gpuErrchk( cudaStreamDestroy(streams[i]) );
    // Free memory
    gpuErrchk( cudaFree(peakBinsDevice[i]) );
    gpuErrchk( cudaFree(peakScoresDevice[i]) );
    gpuErrchk( cudaFree(polyfitsDevice[i]) );
    gpuErrchk( cudaFree(maxScoresDevice[i]) );
    gpuErrchk( cudaFree(pvalsDevice[i]) );
    gpuErrchk( cudaFreeHost(pvalsHost[i]) );
    gpuErrchk( cudaFree(queryPeakBinsDevice[i]) );
  }
}

void synchronizeStream(int streamIdx) {
  cudaStreamSynchronize(streams[streamIdx]);
}

__host__
void runKernel(short *peakBins, short *peakScores, double *polyfits, int *maxScores, size_t N, short *queryPeakBins, size_t M, double** pvalsHost, int streamIdx) {
  cudaSetDevice(0);
  
  cudaStream_t *stream = &streams[streamIdx];
  /*
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  
  cudaEventRecord(start);
  */
  gpuErrchk( cudaMemcpy(peakBinsDevice[streamIdx], peakBins, N * SCORING_PEAKS * sizeof(short), cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(peakScoresDevice[streamIdx], peakScores, N * SCORING_PEAKS * sizeof(short), cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(polyfitsDevice[streamIdx], polyfits, N * POLYFIT_SIZE * sizeof(double), cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(maxScoresDevice[streamIdx], maxScores, N * sizeof(int), cudaMemcpyHostToDevice) );
  gpuErrchk( cudaMemcpy(queryPeakBinsDevice[streamIdx], queryPeakBins, M * SCORING_PEAKS * sizeof(short), cudaMemcpyHostToDevice) );
  
  
  dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
  dim3 dimGrid;
  dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
  dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;  
  
  calculatePvals<<<dimGrid, dimBlock, 0, *stream>>>(peakBinsDevice[streamIdx], peakScoresDevice[streamIdx], polyfitsDevice[streamIdx], maxScoresDevice[streamIdx], N, queryPeakBinsDevice[streamIdx], M, pvalsDevice[streamIdx]);
  
  
  //gpuErrchk( cudaDeviceSynchronize() );
  
  gpuErrchk( cudaMemcpyAsync(pvalsHost[streamIdx], pvalsDevice[streamIdx], N * PVEC_MAX_BATCH_SIZE * sizeof(double), cudaMemcpyDeviceToHost, *stream) );
  /*
  double sumPvals = 0;
  for (size_t i = 0; i < N; i++) {
    sumPvals += pvals[i*M];
  }
  std::cout << "Sum: " << sumPvals << std::endl;
  */
  
  /*
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  
  totaltime += milliseconds;
  
  std::cout << "Elapsed time (ms): " << totaltime << std::endl;
  */
}

} /* namespace maracluster */
