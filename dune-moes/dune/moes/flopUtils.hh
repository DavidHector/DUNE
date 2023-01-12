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

size_t flopsCompGenMinMag(size_t iterations, size_t N, size_t nev, size_t qrFrequency, size_t L, size_t U, size_t M)
{
    size_t gSflops = 2 * nev * nev * N - nev * nev / 2 - nev / 2;
    size_t inverseFlops = 4 * nev * N + 2 * nev * L + 2 * nev * U - nev;
    size_t sparseMatmul = 2 * nev * M - nev * N;
    size_t getEvsflops = 2 * nev * M + 3 * nev * N; // This is not nev and it is also not the same for generalized and standard problems
    size_t gsNumber = iterations / qrFrequency;
    size_t total = gsNumber * gSflops + iterations * inverseFlops + iterations * sparseMatmul + iterations * getEvsflops;
    return total;
}

size_t flopsCompStdMinMag()
{
    //TODO
}

size_t flopsCompGenMaxMag()
{
    //TODO
}

size_t flopsCompStdMaxMag()
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
    size_t iterations, iterationstmp, L, U, Annz;
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
        for (size_t j = 0; j < repetitions[i]; j++)
        {
            moesflops.computeGenMinMagnitudeApprox(B, tolerance, eigenvecs, eigenvals, rhsWidths[i], qrFrequency, sigma, alpha, L, U, LUflops, iterationstmp);
            iterations += iterationstmp;
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        double averageDuration = (double)duration.count() / (double)repetitions[i];
        iterations = iterations / repetitions[i];
        flops = flopsCompGenMinMag(iterations, N, rhsWidths[i], qrFrequency, L, U, Annz) + LUflops;
        gflops = flops / averageDuration;
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
        flops = flopsCompGenMinMag(iterations, N, rhsWidth, qrFrequency, L, U, Annz) + LUflops;
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
    size_t L, U, iterations, Annz;
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
            auto startMoes = std::chrono::high_resolution_clock::now();
            for (size_t reps = 0; reps < repetitions[irhs]; reps++)
            {
                flopsMoes.computeGenMinMagnitude(B, tolerance, moesevs, moeslambdas, rhsWidths[irhs], qrFrequency, sigma, L, U, LUflops, iterations);
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

            flopsM = flopsCompGenMinMag(iterations, Ns[iN], rhsWidths[irhs], qrFrequency, L, U, Annz) + LUflops;
            gflopsM = flopsM / averageDurationM;
            outFile << "\n"
                    << Ns[iN] << "," << rhsWidths[irhs] << "," << repetitions[irhs] << "," << gflopsM << "," << averageDurationM << "," << averageDurationA << ",";
        }
    }

    outFile.close();
}

void flopsParGenMinMag() {}

void flopsSeqGenMaxMag() {}

void flopsParGenMaxMag() {}

// TODO: Compare all sequential things against arpack

#endif