#include <iostream>
#include <vector>
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

int main(int argc, char const *argv[])
{
    size_t N = 10000;
    size_t rhsWidth = 256;
    size_t repetitions = 10;
    const double tolerance = 1e-10;

    if (argc > 1)
    {
        if (1 == std::sscanf(argv[1], "%zu", &N))
        {
        }
        else
        {
            std::cout << "Please enter an unsigned integer!" << std::endl;
            return -1;
        }
    }
    if (argc > 2)
    {
        if (1 == std::sscanf(argv[2], "%zu", &rhsWidth))
        {
        }
        else
        {
            std::cout << "Please enter a power of 32!" << std::endl;
            return -1;
        }
    }
    if (argc > 3)
    {
        if (1 == std::sscanf(argv[3], "%zu", &repetitions))
        {
        }
        else
        {
            std::cout << "Please enter an unsigned integer!" << std::endl;
            return -1;
        }
    }
    size_t Qsize = N * rhsWidth;
    size_t qCols = rhsWidth / 8;
    size_t EVNumber = 8;
    double shift = 0.0;
    std::unique_ptr<double[]> Q(new double[Qsize]);
    std::shared_ptr<double[]> Qs(new double[Qsize]);
    std::vector<double> EVs(EVNumber, 0.0);
    static const int BS = 1;
    typedef Dune::FieldMatrix<double, BS, BS> MatrixBlock;
    typedef Dune::BCRSMatrix<MatrixBlock> BCRSMat; // Matrix Type
    typedef Dune::BlockVector<Dune::FieldVector<double, BS>> VEC;
    BCRSMat laplacian;
    BCRSMat identity;
    BCRSMat identitylap;
    BCRSMat B;
    BCRSMat neumann;
    // setupLaplacian(B, std::sqrt(N));
    setupLaplacian(laplacian, std::sqrt(N)); //AAAAAArgghh, this sets the matrix to size N*N x N*N (with N*N*5 entries, not really sure what the BlockSize does)
    setupIdentity(neumann, N);
    setupIdentity(identity, N);
    std::string slap = "toMatlabWriter.txt";
    Dune::writeMatrixToMatlab(laplacian, slap);
    std::string slapM = "toMatrixMarket.txt";
    Dune::storeMatrixMarket(laplacian, slapM);
    setupIdentityWithLaplacianSparsityPattern(identitylap, N);
    neumann[0][0] = 0.0;
    neumann[N - 1][N - 1] = 0.0;
    transposeMatMultMat(B, neumann, laplacian);
    matMultMat(B, laplacian, neumann);
    ArpackMLGeneo::ArPackPlusPlus_Algorithms<BCRSMat, VEC> arpack(laplacian);
    moes<BCRSMat, VEC> moes(laplacian);
    // How to produce the right hand-side?
    VEC vec(N); //for template
    vec = 0.0;  //just copying from multilevel_geneo_preconditioner.hh
    std::vector<VEC> eigenvecs(EVNumber, vec);
    std::vector<VEC> moeseigenvecs(EVNumber, vec);
    std::vector<double> eigenvals(EVNumber, 0.0);
    std::vector<double> moeseigenvals(EVNumber, 0.0);
    std::cout << "init of arpack was successfull" << std::endl
              << "starting arpack computeGenSymShiftInvertMinMagnitude: " << std::endl;
    // arpack.computeGenSymShiftInvertMinMagnitude(B, tolerance, eigenvecs, eigenvals, shift);
    // can get matrix rows via const int nrows/ncols = A.nrows()/A.ncols()
    // The shift happens inside here with the MAT ashiftb(A) ashiftb.axpy(-sigma, b_) function
    // arpack.computeGenSymShiftInvertMinMagnitudeAdaptive(B,tolerance,eigenvecs,eigenvals,shift,eigenvalue_fine_threshold_,n_eigenvectors_fine_used_);

    // How to get B?
    // moes.computeStdMaxMagnitude(tolerance, moeseigenvecs, eigenvals, EVNumber, 2);
    moes.computeGenMaxMagnitude(identitylap, tolerance, moeseigenvecs, eigenvals, EVNumber, 2, -0.5);
    largestEVsIterative(laplacian, Q, qCols, N, 10000, 1);
    getEigenvalues(laplacian, Q, qCols, N, EVs);
    std::cout << "The largest " << EVNumber << " Eigenvalues are: " << std::endl;
    for (size_t i = 0; i < EVNumber; i++)
    {
        std::cout << EVs[i] << std::endl;
    }

    std::cout << "The largest " << EVNumber << " moes Eigenvalues are: " << std::endl;
    for (size_t i = 0; i < EVNumber; i++)
    {
        std::cout << eigenvals[i] << std::endl;
    }

    // moes.computeStdMinMagnitude(tolerance, moeseigenvecs, moeseigenvals, EVNumber, 2, -0.5); // Gotta use a negative shift, otherwise weird things happen
    size_t L, U, iterations;
    double LUflops;
    moes.computeGenMinMagnitude(identity, tolerance, moeseigenvecs, moeseigenvals, EVNumber, 2, -0.5, L, U, LUflops, iterations);
    // Should do the same as smallest EVs Iterative
    arpack.computeStdNonSymMinMagnitude(identity, tolerance, eigenvecs, eigenvals, -0.5);
    std::vector<VEC> redeigenvecs(4, vec);
    redeigenvecs[0] = eigenvecs[0];
    redeigenvecs[1] = eigenvecs[1];
    redeigenvecs[2] = eigenvecs[2];
    redeigenvecs[3] = eigenvecs[3];
    std::vector<VEC> redmoeseigenvecs(4, vec);
    redmoeseigenvecs[0] = moeseigenvecs[0];
    redmoeseigenvecs[1] = moeseigenvecs[1];
    redmoeseigenvecs[2] = moeseigenvecs[2];
    redmoeseigenvecs[3] = moeseigenvecs[3];
    double csN = moes.columnSumNorm(redmoeseigenvecs, redeigenvecs);

    std::cout << std::endl
              << "The smallest " << EVNumber << " moes Eigenvalues are: " << std::endl;
    for (size_t i = 0; i < EVNumber; i++)
    {
        std::cout << moeseigenvals[i] << std::endl;
    }

    std::cout << std::endl
              << "The smallest " << EVNumber << " arpack Eigenvalues are: " << std::endl;
    for (size_t i = 0; i < EVNumber; i++)
    {
        std::cout << eigenvals[i] << std::endl;
    }
    std::cout << "The column sum norm (for the first four eigenvecs) is: " << csN << std::endl;

    return 0;
}
