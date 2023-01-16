#ifndef DUNE_MOES_FLOPUTILS_HH
#define DUNE_MOES_FLOPUTILS_HH

#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <chrono>
#ifdef HAVE_CONFIG_H // Whelp, this should always be included, but I didnt
#include "config.h"
#endif
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/test/laplacian.hh>
#include <dune/istl/test/identity.hh>
#include <dune/istl/matrixmatrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/common/fmatrix.hh>
#include <dune/moes/MatrixMult.hh>
#include <dune/moes/moes.hh>
#include <dune/moes/Utils.hh>
#include <dune/moes/arpack_geneo_wrapper.hh>

double flopsCompGenMinMagIterationSum(size_t iterations, size_t N, size_t nev, size_t qrFrequency, size_t L, size_t U, size_t M)
{
    double iterationsD = (double)iterations;
    double ND = (double)N;
    double nevD = (double)nev;
    double qrFrequencyD = (double)qrFrequency;
    double LD = (double)L;
    double UD = (double)U;
    double MD = (double)M;

    double gSflops = iterationsD / qrFrequencyD * (2.0 * nevD * nevD * ND - nevD * nevD / 2.0 - nevD / 2.0);
    double inverseFlops = 2.0 * nevD * LD + 2.0 * nevD * UD + 2.0 * nevD * ND - nevD;
    double sparseMatmul = 2.0 * nevD * MD - nevD * ND;
    double getEvsflops = 4.0 * nevD * MD + 2.0 * nevD * ND;
    double total = gSflops + inverseFlops + sparseMatmul + getEvsflops;
    return total;
}

double flopsCompStdMinMag(size_t iterations, size_t N, size_t nev, size_t qrFrequency, size_t L, size_t U, size_t M)
{
    double gSflops = iterations / qrFrequency * (2.0 * nev * nev * N - nev * nev / 2.0 - nev / 2.0);
    double inverseFlops = 2.0 * nev * L + 2.0 * nev * U + 2.0 * nev * N - nev;
    double sparseMatmul = 2.0 * nev * M - nev * N;
    double getEvsflops = 2.0 * nev * M + 3.0 * nev * N;
    double total = gSflops + inverseFlops + sparseMatmul + getEvsflops;
    total *= iterations;
    return total;
}

void flopsCompGenMaxMag()
{
    //TODO
}

void flopsCompStdMaxMag()
{
    //TODO
}

/**
 * @brief Flop measurement for Matrices with kernel intersections (sequential execution)
 * 
 * @tparam MAT 
 * @tparam VEC 
 * @param filenameA 
 * @param filenameB 
 * @param filenameOut 
 * @param tolerance 
 * @param sigma 
 * @param alpha 
 * @param qrFrequency 
 */
template <typename MAT, typename VEC>
void flopsSeqGenMinApproxFileRead(const std::string filenameA, const std::string filenameB, const std::string filenameOut, const double tolerance = 1e-8, const double sigma = 0.01, const double alpha = 0.001, const size_t qrFrequency = 1)
{
    size_t evsout = 8;
    MAT A, B;
    Dune::loadMatrixMarket(A, filenameA);
    Dune::loadMatrixMarket(B, filenameB);
    moes<MAT, VEC> moesflops(A);
    size_t N = A.N();
    VEC vec(N);
    vec = 0.0;
    std::vector<VEC> eigenvecs(evsout, vec);
    std::vector<double> eigenvals(evsout, 0.0);
    const size_t lenRhsWidths = 6;
    size_t rhsWidths[lenRhsWidths] = {8, 16, 24, 32, 40, 48};
    size_t repetitions[lenRhsWidths] = {500, 100, 50, 10, 10, 1};
    size_t iterations, sumIterations, L, U, Annz;
    Annz = A.nonzeroes();
    double gflops, flops;
    double LUflops;

    std::ofstream outputFile;
    outputFile.open(filenameOut);
    outputFile << "rhsWidth, repetitions, iterations, GFLOPs,";

    for (size_t i = 0; i < lenRhsWidths; i++)
    {
        iterations = 0;
        auto start = std::chrono::high_resolution_clock::now();
        sumIterations = 0;
        for (size_t j = 0; j < repetitions[i]; j++)
        {
            moesflops.computeGenMinMagnitudeApprox(B, tolerance, eigenvecs, eigenvals, rhsWidths[i], qrFrequency, sigma, alpha, L, U, LUflops, iterations);
            sumIterations += iterations;
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        double averageDuration = (double)duration.count() / (double)repetitions[i];
        iterations = sumIterations / repetitions[i];
        flops = flopsCompGenMinMagIterationSum(sumIterations, N, rhsWidths[i], qrFrequency, L, U, Annz) + repetitions[i] * LUflops; // all flops I guess
        gflops = flops / duration.count();
        outputFile << "\n"
                   << rhsWidths[i] << "," << repetitions[i] << "," << iterations << "," << gflops << ",";
    }
    outputFile.close();
}

template <typename T>
T vecAvg(std::vector<T> v)
{
    T avg = 0;
    for (auto &&i : v)
    {
        avg += i;
    }
    avg /= v.size();
    return avg;
}

template <typename MAT, typename VEC>
void singleThreadGenMinApprox(MAT &A, MAT &B, const double &epsilon, const int &nev, const int &qrFrequency, const double &sigma, const double &alpha, size_t &L, size_t &U, double &LUflops, size_t &iterations)
{
    moes<MAT, VEC> moesST(A);
    size_t N = A.N();
    VEC vec(N);
    vec = 0.0;
    std::vector<VEC> eigenvecs(8, vec);
    std::vector<double> eigenvals(8, 0.0);
    moesST.computeGenMinMagnitudeApprox(B, epsilon, eigenvecs, eigenvals, nev, qrFrequency, sigma, alpha, L, U, LUflops, iterations);
}

/**
 * @brief Flop measurement for parallely running the genminapprox algorithm when reading from file, solves the problem (A - \sigma B + \alpha I) x = (\lambda - \sigma) B x for smalles evs 
 * 
 * @tparam MAT Matrix Type (Should be BCRSMatrix)
 * @tparam VEC Vector Type
 * @param filenameA File for Matrix A (NOTE: file must be in Matrix Market format, so it needs the header)
 * @param filenameB File for Matrix B (NOTE: file must be in Matrix Market format, so it needs the header)
 * @param filenameOut File to write measurements to (writes in csv format)
 * @param tolerance Tolerance for the solver, i.e. the convergence criterion
 * @param sigma Shifting value for B
 * @param alpha Shifting value for the identity matrix (must be != 0, when trying to approximate the eigenspace for matrices with kernel intersection)
 * @param qrFrequency After how many iterations the rhs Multivector gets orthonormalized
 */
template <typename MAT, typename VEC>
void flopsParGenMinApproxFileRead(const std::string filenameA, const std::string filenameB, const std::string filenameOut, const double tolerance = 1e-8, const double sigma = 0.01, const double alpha = 0.001, const size_t qrFrequency = 1)
{
    size_t evsout = 8;
    MAT A, B;
    Dune::loadMatrixMarket(A, filenameA);
    Dune::loadMatrixMarket(B, filenameB);
    size_t N = A.N();
    size_t rhsWidth = 40;
    const size_t lenThreadCounts = 6;
    size_t threadCounts[lenThreadCounts] = {4, 8, 16, 32, 64, 128};
    size_t repetitions[lenThreadCounts] = {100, 50, 10, 10, 1, 1};
    size_t iterations, iterationstmp, L, U, Annz;
    Annz = A.nonzeroes();
    double gflops, flops;
    double LUflops;

    std::ofstream outputFile;
    outputFile.open(filenameOut);
    outputFile << "rhsWidth, repetitions, iterations, GFLOPs, threadcount,";

    for (size_t i = 0; i < lenThreadCounts; i++)
    {
        iterationstmp = 0.0;
        auto start = std::chrono::high_resolution_clock::now();
        for (size_t j = 0; j < repetitions[i]; j++)
        {
            std::vector<std::thread> threads;
            std::vector<size_t> iterations(threadCounts[i], 0);
            for (size_t tC = 0; tC < threadCounts[i]; tC++)
            {
                threads.push_back(std::thread(singleThreadGenMinApprox<MAT, VEC>, std::ref(A), std::ref(B), std::ref(tolerance), std::ref(rhsWidth), std::ref(qrFrequency), std::ref(sigma), std::ref(alpha), std::ref(L), std::ref(U), std::ref(LUflops), std::ref(iterations[tC])));
            }
            for (size_t tC = 0; tC < threadCounts[i]; tC++)
            {
                threads[tC].join();
            }
            iterationstmp += vecAvg(iterations);
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        double averageDuration = (double)duration.count() / (double)repetitions[i];
        iterations = iterationstmp / repetitions[i];
        flops = flopsCompGenMinMagIterationSum(iterations, N, rhsWidth, qrFrequency, L, U, Annz) + LUflops;
        flops *= threadCounts[i];
        gflops = flops / averageDuration;
        outputFile << "\n"
                   << rhsWidth << "," << repetitions[i] << "," << iterations << "," << gflops << "," << threadCounts[i] << ",";
    }
    outputFile.close();
}

void flopsSeqStdMinMag() {}

void flopsParStdMinMag() {}

void flopsSeqStdMaxMag() {}

void flopsParStdMaxMag() {}

template <typename MAT, typename VEC>
void flopsSeqGenMinMagLap(const std::string filenameOut, const double tolerance = 1e-8, const double sigma = 0.01, const size_t qrFrequency = 5)
{
    const size_t lenNs = 5;
    size_t Ns[5] = {2500, 10000, 40000, 90000, 160000};
    const size_t lenrhsWidths = 7;
    size_t rhsWidths[lenrhsWidths] = {8, 16, 24, 32, 40, 56, 64};
    const size_t lenRepetitions = 7;
    size_t repetitions[lenRepetitions] = {100, 50, 30, 10, 10, 10, 10};
    size_t L, U, iterations, Annz, sumIterations;
    double LUflops, flopsM, gflopsM;

    std::ofstream outFile;
    outFile.open(filenameOut);
    outFile << "N, rhsWidth, repetitions, gflopsmoes[GFLOPS], timeMoes[ns], timeArpack[ns],";
    for (size_t iN = 0; iN < lenNs; iN++)
    {
        std::cout << "N = " << Ns[iN] << std::endl;
        MAT A, B;
        setupLaplacian(A, std::sqrt(Ns[iN]));
        setupIdentity(B, Ns[iN]);
        VEC vec(Ns[iN]);
        vec = 0.0;
        ArpackMLGeneo::ArPackPlusPlus_Algorithms<MAT, VEC> arpack(A);
        moes<MAT, VEC> flopsMoes(A);
        Annz = A.nonzeroes();
        for (size_t irhs = 0; irhs < lenrhsWidths; irhs++)
        {
            std::cout << "rhsWidth = " << rhsWidths[irhs] << std::endl;
            std::vector<VEC> moesevs(rhsWidths[irhs], vec);
            std::vector<VEC> arevs(rhsWidths[irhs], vec);
            std::vector<double> moeslambdas(rhsWidths[irhs], 0.0);
            std::vector<double> arlambdas(rhsWidths[irhs], 0.0);
            sumIterations = 0;
            auto startMoes = std::chrono::high_resolution_clock::now();
            for (size_t reps = 0; reps < repetitions[irhs]; reps++)
            {
                flopsMoes.computeGenMinMagnitude(B, tolerance, moesevs, moeslambdas, rhsWidths[irhs], qrFrequency, sigma, L, U, LUflops, iterations);
                sumIterations += iterations;
            }
            auto stopMoes = std::chrono::high_resolution_clock::now();
            auto durationMoes = std::chrono::duration_cast<std::chrono::nanoseconds>(stopMoes - startMoes);
            double averageDurationM = (double)durationMoes.count() / (double)repetitions[irhs];

            auto startAr = std::chrono::high_resolution_clock::now();
            for (size_t reps = 0; reps < repetitions[irhs]; reps++)
            {
                arpack.computeGenNonSymShiftInvertMinMagnitude(B, tolerance, arevs, arlambdas, sigma);
            }
            auto stopAr = std::chrono::high_resolution_clock::now();
            auto durationAr = std::chrono::duration_cast<std::chrono::nanoseconds>(stopAr - startAr);
            double averageDurationA = (double)durationAr.count() / (double)repetitions[irhs];

            // Possible mistake, I am always using the last iterations number
            flopsM = flopsCompGenMinMagIterationSum(sumIterations, Ns[iN], rhsWidths[irhs], qrFrequency, L, U, Annz) + repetitions[irhs] * LUflops;
            gflopsM = flopsM / durationMoes.count();
            outFile << "\n"
                    << Ns[iN] << "," << rhsWidths[irhs] << "," << repetitions[irhs] << "," << gflopsM << "," << averageDurationM << "," << averageDurationA << ",";
        }
    }

    outFile.close();
}

double flopsGS(const size_t &N, const size_t &W)
{
    return 2.0 * W * W * N - 0.5 * W * W - 0.5 * W;
}

void singleThreadGS(const size_t &N, const size_t &W, const size_t &repetitions)
{
    std::shared_ptr<double[]> Q(new double[N * W]);
    fillMatrixRandom(Q, N * W);
    for (size_t k = 0; k < repetitions; k++)
    {
        qrFixedBlockOptimizedDouble(Q, N, W / 8);
    }
}

void singleThreadGSNaive(const size_t &N, const size_t &W, const size_t &repetitions)
{
    std::shared_ptr<double[]> Q(new double[N * W]);
    fillMatrixRandom(Q, N * W);
    for (size_t k = 0; k < repetitions; k++)
    {
        qrNaiveQNaive(Q, N, W);
    }
}

void flopsGSAutoST(const std::string filenameOut)
{
    const int lenN = 7;
    const int lenrhsWidth = 8;
    size_t Ns[lenN] = {1000, 5000, 10000, 20000, 50000, 100000, 200000};
    size_t repetitions[lenrhsWidth] = {500, 100, 50, 10, 1, 1, 1, 1};
    size_t rhsWidths[lenrhsWidth] = {8, 16, 32, 64, 128, 256, 512, 1024};
    std::ofstream outputFile;
    outputFile.open(filenameOut);
    outputFile << "N, rhsWidth, repetitions, GFLOPS, GFLOPSNaive";
    for (size_t i = 0; i < lenN; i++)
    {
        std::cout << "N = " << Ns[i] << std::endl;
        for (size_t j = 0; j < lenrhsWidth; j++)
        {
            auto start = std::chrono::high_resolution_clock::now();
            singleThreadGS(Ns[i], rhsWidths[j], repetitions[j]);
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            double averageDuration = duration.count() / repetitions[j];
            double gfGS = flopsGS(Ns[i], rhsWidths[j]) / averageDuration;

            start = std::chrono::high_resolution_clock::now();
            singleThreadGSNaive(Ns[i], rhsWidths[j], repetitions[j]);
            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            averageDuration = duration.count() / repetitions[j];
            double gfGSNaive = flopsGS(Ns[i], rhsWidths[j]) / averageDuration;

            outputFile << "\n"
                       << Ns[i] << "," << rhsWidths[j] << "," << repetitions[j] << "," << gfGS << "," << gfGSNaive << ",";
        }
    }
    outputFile.close();
}

void flopsGSAutoMT(const std::string filenameOut)
{
    const int lenN = 7;
    const int lenrhsWidth = 8;
    size_t Ns[lenN] = {1000, 5000, 10000, 20000, 50000, 100000, 200000};
    size_t repetitions[lenrhsWidth] = {500, 100, 50, 10, 1, 1, 1, 1};
    size_t rhsWidths[lenrhsWidth] = {8, 16, 32, 64, 128, 256, 512, 1024};
    const size_t lenThreadCounts = 6;
    size_t threadCounts[lenThreadCounts] = {4, 8, 16, 32, 64, 128};
    std::ofstream outputFile;
    outputFile.open(filenameOut);
    outputFile << "N, rhsWidth, threadCount, repetitions, GFLOPS, GFLOPSNaive,";
    for (size_t i = 0; i < lenN; i++)
    {
        std::cout << "N = " << Ns[i] << std::endl;
        for (size_t j = 0; j < lenrhsWidth; j++)
        {
            for (size_t tC = 0; tC < lenThreadCounts; tC++)
            {
                std::vector<std::thread> threads;
                // Vectorized
                auto start = std::chrono::high_resolution_clock::now();
                for (size_t t = 0; t < threadCounts[tC]; t++)
                {
                    threads.push_back(std::thread(singleThreadGS, std::ref(Ns[i]), std::ref(rhsWidths[j]), std::ref(repetitions[j])));
                }
                for (size_t t = 0; t < threadCounts[tC]; t++)
                {
                    threads[t].join();
                    threads.pop_back();
                }
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
                double averageDuration = duration.count() / (repetitions[j] * threadCounts[tC]);
                double gfGS = flopsGS(Ns[i], rhsWidths[j]) / averageDuration;

                // Non-vectorized
                start = std::chrono::high_resolution_clock::now();
                for (size_t t = 0; t < threadCounts[tC]; t++)
                {
                    threads.push_back(std::thread(singleThreadGSNaive, std::ref(Ns[i]), std::ref(rhsWidths[j]), std::ref(repetitions[j])));
                }
                for (size_t t = 0; t < threadCounts[tC]; t++)
                {
                    threads[t].join();
                    threads.pop_back();
                }
                stop = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
                averageDuration = duration.count() / (repetitions[j] * threadCounts[tC]);
                double gfGSNaive = flopsGS(Ns[i], rhsWidths[j]) / averageDuration;

                outputFile << "\n"
                           << Ns[i] << "," << rhsWidths[j] << "," << threadCounts[tC] << "," << repetitions[j] << "," << gfGS << "," << gfGSNaive << ",";
            }
        }
    }
    outputFile.close();
}

template <typename MAT>
double flopsMatmul(const MAT &A, const size_t &W)
{
    size_t N = A.N();
    size_t nnz = A.nonzeroes();
    return 2 * W * nnz - W * N;
}

template <typename MAT>
void singleThreadMatmul(const MAT &A, const size_t &W, const size_t &repetitions)
{
    const size_t N = A.N();
    std::shared_ptr<double[]> Q(new double[N * W]);
    std::shared_ptr<double[]> AQ(new double[N * W]);
    fillMatrixRandom(Q, N * W);
    for (size_t k = 0; k < repetitions; k++)
    {
        MultQSimple(A, Q, AQ, W / 8, N);
    }
}

template <typename MAT>
void singleThreadMatmulNaive(const MAT &A, const size_t &W, const size_t &repetitions)
{
    const size_t N = A.N();
    std::shared_ptr<double[]> Q(new double[N * W]);
    std::shared_ptr<double[]> AQ(new double[N * W]);
    fillMatrixRandom(Q, N * W);
    for (size_t k = 0; k < repetitions; k++)
    {
        MultQSimpleNaiveQNaive<MAT>(A, Q, AQ, W, N);
    }
}

void flopsMatmulST(const std::string filenameOut)
{
    const int BS = 1;
    typedef Dune::FieldMatrix<double, BS, BS> MatrixBlock;
    typedef Dune::BCRSMatrix<MatrixBlock> BCRSMat;
    const int lenN = 6;
    const int lenrhsWidth = 8;
    size_t Ns[lenN] = {10000, 40000, 90000, 160000, 490000, 1000000};
    size_t repetitions[lenN] = {50, 10, 5, 1, 1};
    size_t rhsWidths[lenrhsWidth] = {8, 16, 32, 64, 128, 256, 512, 1024};
    std::ofstream outputFile;
    outputFile.open(filenameOut);
    outputFile << "N, rhsWidth, repetitions, GFIdentity, GFIdentityNaive, GFLaplacian, GFLaplacianNaive";
    for (size_t i = 0; i < lenN; i++)
    {
        std::cout << "N = " << Ns[i] << std::endl;
        for (size_t j = 0; j < lenrhsWidth; j++)
        {
            BCRSMat identity;
            BCRSMat laplacian;
            setupIdentity(identity, Ns[i]);
            setupLaplacian(laplacian, std::sqrt(Ns[i]));

            // matrix multiplication identity matrix, vectorized
            auto start = std::chrono::high_resolution_clock::now();
            singleThreadMatmul(identity, rhsWidths[j], repetitions[i]);
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            double averageDuration = duration.count() / repetitions[i];
            double gfIdentity = flopsMatmul(identity, rhsWidths[j]) / averageDuration;

            // matrix multiplication identity matrix, naive
            start = std::chrono::high_resolution_clock::now();
            singleThreadMatmulNaive(identity, rhsWidths[j], repetitions[i]);
            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            averageDuration = duration.count() / repetitions[i];
            double gfIdentityNaive = flopsMatmul(identity, rhsWidths[j]) / averageDuration;

            // matrix multiplication laplacian matrix, vectorized
            start = std::chrono::high_resolution_clock::now();
            singleThreadMatmul(laplacian, rhsWidths[j], repetitions[i]);
            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            averageDuration = duration.count() / repetitions[i];
            double gfLaplacian = flopsMatmul(laplacian, rhsWidths[j]) / averageDuration;

            // matrix multiplication laplacian matrix, naive
            start = std::chrono::high_resolution_clock::now();
            singleThreadMatmulNaive(identity, rhsWidths[j], repetitions[i]);
            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            averageDuration = duration.count() / repetitions[i];
            double gfLaplacianNaive = flopsMatmul(laplacian, rhsWidths[j]) / averageDuration;

            outputFile << "\n"
                       << Ns[i] << "," << rhsWidths[j] << "," << repetitions[i] << "," << gfIdentity << "," << gfIdentityNaive << "," << gfLaplacian << "," << gfLaplacianNaive << ",";
        }
    }
    outputFile.close();
}

void flopsMatmulMT(const std::string filenameOut)
{
    const int BS = 1;
    typedef Dune::FieldMatrix<double, BS, BS> MatrixBlock;
    typedef Dune::BCRSMatrix<MatrixBlock> BCRSMat;
    const int lenN = 5;
    const int lenrhsWidth = 8;
    size_t Ns[lenN] = {10000, 40000, 160000, 490000, 1000000};
    size_t repetitions[lenN] = {50, 10, 5, 1, 1};
    size_t rhsWidths[lenrhsWidth] = {8, 16, 32, 64, 128, 256, 512, 1024};
    const size_t lenThreadCounts = 6;
    size_t threadCounts[lenThreadCounts] = {4, 8, 16, 32, 64, 128};

    std::ofstream outputFile;
    outputFile.open(filenameOut);
    outputFile << "N, rhsWidth, repetitions, threadCount, GFIdentity, GFIdentityNaive, GFLaplacian, GFLaplacianNaive";
    for (size_t i = 0; i < lenN; i++)
    {
        std::cout << "N = " << Ns[i] << std::endl;
        for (size_t j = 0; j < lenrhsWidth; j++)
        {
            BCRSMat identity;
            BCRSMat laplacian;
            setupIdentity(identity, Ns[i]);
            setupLaplacian(laplacian, std::sqrt(Ns[i]));

            for (size_t tC = 0; tC < lenThreadCounts; tC++)
            {
                std::vector<std::thread> threads;
                // matrix multiplication identity matrix, vectorized
                auto start = std::chrono::high_resolution_clock::now();
                for (size_t t = 0; t < threadCounts[tC]; t++)
                {
                    threads.push_back(std::thread(singleThreadMatmul<BCRSMat>, std::ref(identity), std::ref(rhsWidths[j]), std::ref(repetitions[i])));
                }
                for (size_t t = 0; t < threadCounts[tC]; t++)
                {
                    threads[t].join();
                    threads.pop_back();
                }
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
                double averageDuration = duration.count() / (repetitions[i] * threadCounts[tC]);
                double gfIdentity = flopsMatmul(identity, rhsWidths[j]) / averageDuration;

                // matrix multiplication identity matrix, naive
                start = std::chrono::high_resolution_clock::now();
                for (size_t t = 0; t < threadCounts[tC]; t++)
                {
                    threads.push_back(std::thread(singleThreadMatmulNaive<BCRSMat>, std::ref(identity), std::ref(rhsWidths[j]), std::ref(repetitions[i])));
                }
                for (size_t t = 0; t < threadCounts[tC]; t++)
                {
                    threads[t].join();
                    threads.pop_back();
                }
                stop = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
                averageDuration = duration.count() / (repetitions[i] * threadCounts[tC]);
                double gfIdentityNaive = flopsMatmul(identity, rhsWidths[j]) / averageDuration;

                // matrix multiplication laplacian matrix, vectorized
                start = std::chrono::high_resolution_clock::now();
                for (size_t t = 0; t < threadCounts[tC]; t++)
                {
                    threads.push_back(std::thread(singleThreadMatmul<BCRSMat>, std::ref(laplacian), std::ref(rhsWidths[j]), std::ref(repetitions[i])));
                }
                for (size_t t = 0; t < threadCounts[tC]; t++)
                {
                    threads[t].join();
                    threads.pop_back();
                }
                stop = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
                averageDuration = duration.count() / (repetitions[i] * threadCounts[tC]);
                double gfLaplacian = flopsMatmul(laplacian, rhsWidths[j]) / averageDuration;

                // matrix multiplication laplacian matrix, naive
                start = std::chrono::high_resolution_clock::now();
                for (size_t t = 0; t < threadCounts[tC]; t++)
                {
                    threads.push_back(std::thread(singleThreadMatmulNaive<BCRSMat>, std::ref(laplacian), std::ref(rhsWidths[j]), std::ref(repetitions[i])));
                }
                for (size_t t = 0; t < threadCounts[tC]; t++)
                {
                    threads[t].join();
                    threads.pop_back();
                }
                stop = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
                averageDuration = duration.count() / (repetitions[i] * threadCounts[tC]);
                double gfLaplacianNaive = flopsMatmul(laplacian, rhsWidths[j]) / averageDuration;

                outputFile << "\n"
                           << Ns[i] << "," << rhsWidths[j] << "," << repetitions[i] << "," << threadCounts[tC] << "," << gfIdentity << "," << gfIdentityNaive << "," << gfLaplacian << "," << gfLaplacianNaive << ",";
            }
        }
    }
    outputFile.close();
}

void flopsParGenMinMag() {}

void flopsSeqGenMaxMag() {}

void flopsParGenMaxMag() {}

// TODO: Compare all sequential things against arpack

#endif