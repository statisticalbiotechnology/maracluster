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
void binaryMatchPeakBinsOld(short *targetPeakBins, short *targetPeakScores, short *queryPeakBins, int *score) {
  size_t qIdx = 0;
  for (size_t i = 0; i < SCORING_PEAKS; ++i) {
    if (qIdx >= SCORING_PEAKS) break;
    while (targetPeakBins[i] > queryPeakBins[qIdx]) {
      if (++qIdx >= SCORING_PEAKS) break;
    }
    if (qIdx >= SCORING_PEAKS) break;
    if (targetPeakBins[i] == queryPeakBins[qIdx]) {
      ++qIdx;
      *score += targetPeakScores[i];
    }
  }
}

__device__
void binaryMatchPeakBins(short *targetPeakBins, short *targetPeakScores, short *queryPeakBins, int *score) {
  size_t tIdx = 0, qIdx = 0;
  short t = -1, q = -1;
  while (tIdx < SCORING_PEAKS && qIdx < SCORING_PEAKS) {
    t = targetPeakBins[tIdx];
    q = queryPeakBins[qIdx];
    if (t == q) {
      *score += targetPeakScores[tIdx];
      ++tIdx;
      ++qIdx;
    } else {
      tIdx += (t < q);
      qIdx += (q < t);
    }
  }
}

__global__
void calculatePvals(short *peakBins, short *peakScores, double *polyfits, int *maxScores, int n, short *queryPeakBins, int m, double *pvals) {
  __shared__ short peakBinsShared[BLOCK_SIZE][SCORING_PEAKS];
  __shared__ short peakScoresShared[BLOCK_SIZE][SCORING_PEAKS];
  __shared__ double polyfitsShared[BLOCK_SIZE][POLYFIT_SIZE];
  __shared__ int maxScoresShared[BLOCK_SIZE];
  __shared__ short queryPeakBinsShared[BLOCK_SIZE][SCORING_PEAKS];
  
  size_t targetIdx = blockIdx.x * BLOCK_SIZE + threadIdx.x;
  size_t queryIdx = blockIdx.y * BLOCK_SIZE + threadIdx.y;
  
  size_t flattenedTargetIdx = targetIdx * SCORING_PEAKS;
  size_t flattenedQueryIdx = queryIdx * SCORING_PEAKS;
  
  size_t placeHolderQueryIdx = blockIdx.y * BLOCK_SIZE + threadIdx.x;
  size_t placeholderFlattenedQueryIdx = placeHolderQueryIdx * SCORING_PEAKS;
  
  if (threadIdx.y == 0 && targetIdx < n) {
    for (int i = 0; i < SCORING_PEAKS; ++i) {
      peakBinsShared[threadIdx.x][i] = peakBins[flattenedTargetIdx + i];
    }
  }
  if (threadIdx.y == 2 && targetIdx < n) {
    for (int i = 0; i < SCORING_PEAKS; ++i) {
      peakScoresShared[threadIdx.x][i] = peakScores[flattenedTargetIdx + i];
    }
  }
  if (threadIdx.y == 4 && targetIdx < n) {
    for (int i = 0; i < POLYFIT_SIZE; ++i) {
      polyfitsShared[threadIdx.x][i] = polyfits[targetIdx * POLYFIT_SIZE + i];
    }
    maxScoresShared[threadIdx.x] = maxScores[targetIdx];
  }
  if (threadIdx.y == 6 && placeHolderQueryIdx < m) {
    for (int i = 0; i < SCORING_PEAKS; ++i) {
      queryPeakBinsShared[threadIdx.x][i] = queryPeakBins[placeholderFlattenedQueryIdx + i];
    }
  }
  __syncthreads();
  if (targetIdx < n && queryIdx < m) {  
    int score = 0;
    binaryMatchPeakBins(&peakBinsShared[threadIdx.x][0], &peakScoresShared[threadIdx.x][0], &queryPeakBinsShared[threadIdx.y][0], &score);
    //binaryMatchPeakBins(&peakBinsShared[threadIdx.x][0], &peakScores[flattenedTargetIdx], &queryPeakBinsShared[threadIdx.y][0], &score);
    //binaryMatchPeakBins(&peakBinsShared[threadIdx.x][0], &peakScores[flattenedTargetIdx], &queryPeakBins[flattenedQueryIdx], &score);
    //binaryMatchPeakBins(&peakBins[flattenedTargetIdx], &peakScores[flattenedTargetIdx], &queryPeakBins[flattenedQueryIdx], &score);
    
    double relScore = static_cast<double>(maxScoresShared[threadIdx.x] - score) / maxScoresShared[threadIdx.x];
    //double relScore = static_cast<double>(maxScores[targetIdx] - score) / maxScores[targetIdx];
    //double relScore = static_cast<double>(score) / maxScores[targetIdx];
    polyval(&polyfitsShared[threadIdx.x][0], &relScore);
    //polyval(&polyfits[targetIdx * POLYFIT_SIZE], &relScore);
    size_t pvalIdx = targetIdx * PVEC_MAX_BATCH_SIZE + queryIdx;
    pvals[pvalIdx] = static_cast<double>(relScore);
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
  int nDevices;

  cudaGetDeviceCount(&nDevices);
  std::cerr << "#GPUs: " << nDevices << std::endl;
  for (int i = 0; i < NUM_STREAMS; ++i) { 
    gpuErrchk( cudaSetDevice(i % NUM_DEVICES) );
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
    gpuErrchk( cudaSetDevice(i % NUM_DEVICES) );
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
  gpuErrchk( cudaSetDevice(streamIdx % NUM_DEVICES) );
  
  cudaStream_t *stream = &streams[streamIdx];
  /*
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  */
  
  /* These cudaMemcpyAsync probably do not work exactly as expected, as the host memory is not pinned. 
     However, it needs to be on the same stream as the kernel to ensure the copies are finished before the kernel executes, which occurs when using cudaMemcpy with data transfers below 64KB */
  gpuErrchk( cudaMemcpyAsync(peakBinsDevice[streamIdx], peakBins, N * SCORING_PEAKS * sizeof(short), cudaMemcpyHostToDevice, *stream) );
  gpuErrchk( cudaMemcpyAsync(peakScoresDevice[streamIdx], peakScores, N * SCORING_PEAKS * sizeof(short), cudaMemcpyHostToDevice, *stream) );
  gpuErrchk( cudaMemcpyAsync(polyfitsDevice[streamIdx], polyfits, N * POLYFIT_SIZE * sizeof(double), cudaMemcpyHostToDevice, *stream) );
  gpuErrchk( cudaMemcpyAsync(maxScoresDevice[streamIdx], maxScores, N * sizeof(int), cudaMemcpyHostToDevice, *stream) );
  gpuErrchk( cudaMemcpyAsync(queryPeakBinsDevice[streamIdx], queryPeakBins, M * SCORING_PEAKS * sizeof(short), cudaMemcpyHostToDevice, *stream) );
  
  //gpuErrchk( cudaDeviceSynchronize() );
  
  dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
  dim3 dimGrid;
  dimGrid.x = (N + dimBlock.x - 1) / dimBlock.x;
  dimGrid.y = (M + dimBlock.y - 1) / dimBlock.y;  
  
  
  calculatePvals<<<dimGrid, dimBlock, 0, *stream>>>(peakBinsDevice[streamIdx], peakScoresDevice[streamIdx], polyfitsDevice[streamIdx], maxScoresDevice[streamIdx], N, queryPeakBinsDevice[streamIdx], M, pvalsDevice[streamIdx]); 
  
  
  //gpuErrchk( cudaDeviceSynchronize() );
  //cudaEventRecord(start);
  gpuErrchk( cudaMemcpyAsync(pvalsHost[streamIdx], pvalsDevice[streamIdx], N * PVEC_MAX_BATCH_SIZE * sizeof(double), cudaMemcpyDeviceToHost, *stream) );
  /*
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  
  totaltime += milliseconds;
  
  std::cout << "Elapsed time (ms): " << totaltime << std::endl;
  */
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
