#include <iostream>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/test/laplacian.hh>
#include <dune/moes/MatrixMult.hh>

int main(int argc, char const *argv[])
{
    size_t N = 1000;
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
    double *Q = new double[Qsize];
    static const int BS = 4;
    typedef Dune::FieldMatrix<double, BS, BS> MatrixBlock;
    typedef Dune::BCRSMatrix<MatrixBlock> BCRSMat;
    BCRSMat laplacian;
    setupLaplacian(laplacian, N / BS);

    return 0;
}
