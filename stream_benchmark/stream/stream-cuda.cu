#include <chrono>

#include "../util.h"
#include "stream-util.h"

__global__ void stream(size_t nx, const double *__restrict__ src, double *__restrict__ dest) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i<nx)
        dest[i] = src[i] + 1;
}

int main(int argc, char *argv[]) {
    size_t nx, nItWarmUp, nIt;
    parseCLA_1d(argc, argv, nx, nItWarmUp, nIt);

    //allocate memory

    size_t size = sizeof(double) * nx;

    double *src, *dest;
    cudaMallocHost(&src, size);
    cudaMallocHost(&dest, size);
    
    double *d_src, *d_dest;
    cudaMalloc(&d_src, size);
    cudaMalloc(&d_dest, size);

    
    // init
    initStream(src, nx);
    // copy from cpu to gpu
    cudaMemcpy(d_src, src, size, cudaMemcpyHostToDevice);

    auto numThreadsPerBlock = 64 ;
    auto numBlocks = (nx+ numThreadsPerBlock-1) / numThreadsPerBlock;

    // warm-up
    for (int i = 0; i < nItWarmUp; ++i) {
        stream<<<numBlocks, numThreadsPerBlock>>>(nx, d_src, d_dest);
        cudaDeviceSynchronize();
        std::swap(d_src, d_dest);
    }

    // measurement
    auto start = std::chrono::steady_clock::now();

    for (int i = 0; i < nIt; ++i) {
        stream<<<numBlocks, numThreadsPerBlock>>>(nx, d_src, d_dest);
        cudaDeviceSynchronize();
        std::swap(d_src, d_dest);
    }

    cudaDeviceSynchronize();
    auto end = std::chrono::steady_clock::now();

    //copy from GPU to CPU --> we swapped d_dest with d_src.. thats why we use d_src here
    cudaMemcpy(dest, d_src, size, cudaMemcpyDeviceToHost);

    printStats(end - start, nx, nIt, streamNumReads, streamNumWrites);

    // check solution --> the result is in dest
    checkSolutionStream(dest, nx, nIt + nItWarmUp);

    cudaFree(d_src);
    cudaFree(d_dest);

    cudaFreeHost(src);
    cudaFreeHost(dest);

    return 0;
}
