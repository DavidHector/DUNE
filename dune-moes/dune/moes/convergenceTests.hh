#ifndef DUNE_MOES_CONVERGENCETESTS_HH
#define DUNE_MOES_CONVERGENCETESTS_HH

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

template <typename MAT, typename VEC>
void csnGenMinApprox(const std::string Afile, const std::string Bfile, const std::string filenameOut, const double tolerance = 1e-8, const double sigma = 0.01, const double alpha = 0.001, const size_t qrFrequency = 1)
{
    MAT A, B;
    Dune::loadMatrixMarket(A, Afile);
    Dune::loadMatrixMarket(B, Bfile);
    const size_t lenIterations = 11;
    const size_t rhsWidth = 40;
    size_t iterations[lenIterations] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    VEC vec(A.N());
    vec = 0.0;
    std::vector<VEC> evs(rhsWidth, vec);
    std::vector<VEC> arpacksevs(rhsWidth, vec);
    std::vector<double> lambdas(rhsWidth, 0.0);
    std::vector<double> arlambdas(rhsWidth, 0.0);
    double LUflops, csn;
    size_t L, U;
    moes<MAT, VEC> csnMoes(A);
    MAT identity;
    MAT Aslightshift(A);
    setupIdentity(identity, A.N());
    Aslightshift.axpy(alpha, identity);
    ArpackMLGeneo::ArPackPlusPlus_Algorithms<MAT, VEC> arpack(Aslightshift);
    arpack.computeGenNonSymShiftInvertMinMagnitude(B, tolerance, arpacksevs, arlambdas, sigma);
    std::ofstream outputFile;
    outputFile.open(filenameOut);
    outputFile << "rhsWidth, iterations, columnSumNorm, arpackTolerance, ";
    for (size_t i = 0; i < lenIterations; i++)
    {
        csnMoes.computeGenMinMagnitudeApproxIterations(B, evs, lambdas, rhsWidth * 2, qrFrequency, sigma, alpha, L, U, LUflops, iterations[i]);
        csn = csnMoes.columnSumNorm(evs, arpacksevs);
        outputFile << "\n"
                   << rhsWidth << "," << iterations[i] << "," << csn << "," << tolerance << ",";
    }
    outputFile.close();
}

template <typename VEC>
void reducevectorlength(const std::vector<VEC> &vold, std::vector<VEC> &vnew)
{
    for (size_t i = 0; i < vnew.size(); i++)
    {
        vnew[i] = vold[i];
    }
}

template <typename MAT, typename VEC>
void csnGenMinLap(const std::string filenameOut, const double tolerance = 1e-8, const double sigma = 0.01, const size_t qrFrequency = 1)
{
    const size_t lenNs = 5;
    size_t Ns[lenNs] = {144, 400, 900, 1600, 2500}; // Must be square numbers
    const size_t lenIterations = 12;
    size_t iterations[lenIterations] = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096};
    const size_t lenrhsWidths = 4;
    size_t rhsWidths[lenrhsWidths] = {8, 16, 32, 64}; // Must be multiples of 8
    double csn;
    size_t nonred;
    std::ofstream outputFile;
    outputFile.open(filenameOut);
    outputFile << "N, rhsWidth, iterations, columnSumNorm, arpackTolerance, ";
    for (size_t iN = 0; iN < lenNs; iN++)
    {
        MAT A, B;
        setupLaplacian(A, std::sqrt(Ns[iN]));
        setupIdentity(B, Ns[iN]);
        VEC vec(Ns[iN]);
        vec = 0.0;
        ArpackMLGeneo::ArPackPlusPlus_Algorithms<MAT, VEC> arpack(A);
        moes<MAT, VEC> csnMoes(A);
        for (size_t irhsW = 0; irhsW < lenrhsWidths; irhsW++)
        {
            nonred = rhsWidths[irhsW] * 2;
            std::vector<VEC> moesevs(rhsWidths[irhsW], vec);
            std::vector<VEC> arevs(nonred, vec);
            std::vector<VEC> arevsred(rhsWidths[irhsW], vec);
            std::vector<double> moeslambdas(rhsWidths[irhsW], 0.0);
            std::vector<double> arlambdas(nonred, 0.0);
            arpack.computeGenNonSymShiftInvertMinMagnitude(B, tolerance, arevs, arlambdas, sigma);
            for (size_t iIts = 0; iIts < lenIterations; iIts++)
            {
                csnMoes.computeGenMinMagnitudeIterations(B, moesevs, moeslambdas, iterations[iIts], nonred, qrFrequency, sigma);
                // could also do flop calculation here
                reducevectorlength(arevs, arevsred);
                csn = csnMoes.columnSumNorm(moesevs, arevsred);
                outputFile << "\n"
                           << Ns[iN] << "," << rhsWidths[irhsW] << "," << iterations[iIts] << "," << csn << "," << tolerance << ",";
            }
        }
        std::cout << "N = " << Ns[iN] << " . Finished." << std::endl;
    }
    outputFile.close();
}
#endif // convergenceTests