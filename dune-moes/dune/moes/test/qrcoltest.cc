#include <iostream>
#include <chrono>

#include <dune/moes/qrcol.hh>
#include <dune/moes/vectorclass/vectorclass.h>

double flopsQR(const size_t& N, const size_t& W){
    return W*(4.0*N) + 0.5 * W*(W-1) * (2.0*N + 1.0 + 2.0*N);
}

double memBandwidthGB(const size_t&N, const size_t&W, const double averageTimeNs){
    double memoryReq = N*W*8.0;
    double bW = memoryReq/averageTimeNs; // Bytes/nanosecond = Gigabytes/second
    return bW;
}

template<typename QR, typename F>
void getGFLOPSqr(Vec4d* Q, QR qr, F f, size_t N, size_t rhsWidth, size_t repetitions, size_t blockSize = 2, size_t uBlockSize = 1){
    double flops = f(N, rhsWidth);
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < repetitions; i++)
    {
        qr(Q, N, rhsWidth/blockSize/4, blockSize, uBlockSize);
    }
    std::cout << "Algorithm finished! Calculating metrics..." << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    auto averageDuration = (double) duration.count()/repetitions;
    auto gFlops = (flops/averageDuration) / 1.0; // 1e6 when duration in ms, 1e3 when using µs
    std::cout << "Average Duration: " << averageDuration/1e6 << "ms" << std::endl;
    std::cout << "Memory usage (full Matrix): " << N*rhsWidth*8.0/1e6 << "MB" << std::endl; 
    std::cout << "Memory Bandwidth: " << memBandwidthGB(N, rhsWidth, averageDuration) << "GB/s (max: 25.6GB/s)" << std::endl;
    std::cout << "GFLOPS: " << gFlops << std::endl;
}

template<typename QR, typename F>
void getGFLOPSqr(double* Q, QR qr, F f, size_t N, size_t rhsWidth, size_t repetitions, size_t blockSize = 2, size_t uBlockSize = 1){
    double flops = f(N, rhsWidth);
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < repetitions; i++)
    {
        qr(Q, N, rhsWidth/blockSize/4, blockSize, uBlockSize);
    }
    std::cout << "Algorithm finished! Calculating metrics..." << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    auto averageDuration = (double) duration.count()/repetitions;
    auto gFlops = (flops/averageDuration) / 1.0; // 1e6 when duration in ms, 1e3 when using µs
    std::cout << "Average Duration: " << averageDuration/1e6 << "ms" << std::endl;
    std::cout << "Memory usage (full Matrix): " << N*rhsWidth*8.0/1e6 << "MB" << std::endl; 
    std::cout << "Memory Bandwidth: " << memBandwidthGB(N, rhsWidth, averageDuration) << "GB/s (max: 25.6GB/s)" << std::endl;
    std::cout << "GFLOPS: " << gFlops << std::endl;
}

int main(int argc, char const *argv[])
{
    size_t N = 1000;
    size_t rhsWidth = 256;
    size_t repetitions = 10;
    size_t blockSize = 2;
    const double tolerance = 1e-6;

    if (argc > 1)
    {
        if (1 == std::sscanf(argv[1], "%zu", &N))
        {
        } else {
            std::cout << "Please enter an unsigned integer!" << std::endl;
            return -1;
        }
    }
    if (argc > 2)
    {
        if (1 == std::sscanf(argv[2], "%zu", &rhsWidth))
        {
        } else {
            std::cout << "Please enter a power of 32!" << std::endl;
            return -1;
        }
    }
    if (argc > 3)
    {
        if (1 == std::sscanf(argv[3], "%zu", &repetitions))
        {
        } else {
            std::cout << "Please enter an unsigned integer!" << std::endl;
            return -1;
        }
    }
    size_t matrixSize = N*rhsWidth/4;
    size_t matrixSizeDouble = N*rhsWidth;
    Vec4d* Q = new Vec4d[matrixSize]; // The matrix
    double* Qdouble = new double[matrixSizeDouble];
    std::cout << "Matrix Size: " << matrixSize << std::endl;

    /* 
    std::cout << "QR-algorithm (Vectorized, variable block Size): " << std::endl;
    fillMatrixRandom(Q, matrixSize);
    getGFLOPSqr(Q, qr, flopsQR, N, rhsWidth, repetitions);
    checkOrthoNormality(Q, N, rhsWidth/blockSize/4, blockSize, tolerance);
    std::cout << std::endl;
    */

    std::cout << "QR-algorithm (Vectorized, fixed block Size): " << std::endl;
    fillMatrixRandom(Q, matrixSize);
    getGFLOPSqr(Q, qrunblocked, flopsQR, N, rhsWidth, repetitions);
    checkOrthoNormalityFixed(Q, N, rhsWidth/8, tolerance);
    std::cout << std::endl;

    std::cout << "QR-algorithm (Vectorized, fixed block Size, Optimized): " << std::endl;
    fillMatrixRandom(Q, matrixSize);
    getGFLOPSqr(Q, qrFixedBlockOptimizedVec4d, flopsQR, N, rhsWidth, repetitions);
    checkOrthoNormalityFixed(Q, N, rhsWidth/8, tolerance);
    std::cout << std::endl;

    std::cout << "QR-algorithm (Vectorized, fixed block Size, Optimized, Doubles): " << std::endl;
    fillMatrixRandom(Qdouble, matrixSizeDouble);
    getGFLOPSqr(Qdouble, qrFixedBlockOptimizedDouble, flopsQR, N, rhsWidth, repetitions);
    checkOrthoNormalityFixed(Q, N, rhsWidth/8, tolerance);
    std::cout << std::endl;

    return 0;
}