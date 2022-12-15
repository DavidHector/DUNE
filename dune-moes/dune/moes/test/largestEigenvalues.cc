#include <iostream>
#include <vector>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/test/laplacian.hh>
#include <dune/moes/MatrixMult.hh>
#include <dune/moes/moes.hh>

int main(int argc, char const *argv[])
{
    size_t N = 10000;
    size_t rhsWidth = 256;
    size_t repetitions = 10;
    const double tolerance = 1e-6;

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
    std::unique_ptr<double[]> Q(new double[Qsize]);
    std::shared_ptr<double[]> Qs(new double[Qsize]);
    std::vector<double> EVs(EVNumber, 0.0);
    static const int BS = 1;
    typedef Dune::FieldMatrix<double, BS, BS> MatrixBlock;
    typedef Dune::BCRSMatrix<MatrixBlock> BCRSMat;
    BCRSMat laplacian;
    setupLaplacian(laplacian, std::sqrt(N)); //AAAAAArgghh, this sets the matrix to size N*N x N*N (with N*N*5 entries, not really sure what the BlockSize does)
    largestEVsIterative(laplacian, Q, qCols, N, 1000, 1);
    // std::cout << "After calling largestEVs: Q[0] = " << Q[0] << std::endl;
    getEigenvalues(laplacian, Q, qCols, N, EVs);
    // std::cout << "After calling getEigenvalues: Q[0] = " << Q[0] << std::endl;
    std::cout << "The largest " << EVNumber << " Eigenvalues are: " << std::endl;
    for (size_t i = 0; i < EVNumber; i++)
    {
        std::cout << EVs[i] << std::endl;
    }
    // std::cout << "Before calling smallestEVs: Qs[0] = " << Qs[0] << std::endl;
    // Why does it produce nans after a while?
    smallestEVsIterative(laplacian, Qs, qCols, N, 200, 1);
    // std::cout << "After calling smallestEVs: Qs[0] = " << Qs[0] << std::endl;
    getEigenvalues(laplacian, Qs, qCols, N, EVs);
    // std::cout << "After calling getEigenvalues: Qs[0] = " << Qs[0] << std::endl;
    std::cout << std::endl
              << "The smallest " << EVNumber << " Eigenvalues are: " << std::endl;
    for (size_t i = 0; i < EVNumber; i++)
    {
        std::cout << EVs[i] << std::endl;
    }

    return 0;
}
