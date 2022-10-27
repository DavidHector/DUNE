#include <chrono>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/test/identity.hh>
#include <dune/moes/MatrixMult.hh>
#include <dune/moes/qrcol.hh>

void checkEquality(Vec4d *Qold, Vec4d *Qnew, int matrixSize)
{
    Vec4ib equal;
    for (size_t i = 0; i < matrixSize; i++)
    {
        if (!horizontal_and(Qold[i] == Qnew[i]))
        {
            std::cout << "Vectors are different! Index = " << i << std::endl;
            return;
        }
    }
}

void checkEquality(const double *Qold, const double *Qnew, const int matrixSize)
{
    Vec4ib equal;
    for (size_t i = 0; i < matrixSize; i++)
    {
        if (Qold[i] != Qnew[i])
        {
            std::cout << "Vectors are different! Index = " << i << std::endl;
            return;
        }
    }
}

void printQ(Vec4d *Q, int matrixSize)
{
    std::cout << "Printing Q: " << std::endl;
    for (size_t i = 0; i < matrixSize; i++)
    {
        std::cout << "[" << Q[i][0] << "," << Q[i][1] << "," << Q[i][2] << "," << Q[i][3] << "]"
                  << "\t";
    }
    std::cout << std::endl;
}

int main(int argc, char const *argv[])
{
    // Make a Block Matrix
    int N = 1e4; //
    double Nd = 1e4;
    int rhsWidth = 256;
    int rhsWidthD = 256;
    int matrixSize = N * rhsWidth; // / 4;
    size_t repetitions = 1000;
    double flops = Nd * 4 * rhsWidthD * 2.0 * repetitions; // N matrix blocks, 16 entries per block, 2 operations, rhsWidth vectors
    static const int BS = 4;
    typedef Dune::FieldMatrix<double, BS, BS> MatrixBlock;
    typedef Dune::BCRSMatrix<MatrixBlock> BCRSMat;
    BCRSMat identity;
    setupIdentity(identity, N / BS);

    // Make Qold and Qnew
    /*
    Vec4d *Qold = new Vec4d[matrixSize];
    Vec4d *Qnew = new Vec4d[matrixSize];
    */
    double *Qold = new double[matrixSize];
    double *Qnew = new double[matrixSize];

    fillMatrixRandom(Qold, matrixSize);
    // std::cout << Qold[0][0] << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < repetitions; i++)
    {
        MultQ(identity, Qold, Qnew, rhsWidth / 8, N);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    // std::cout << Qnew[0][0] << std::endl;
    checkEquality(Qold, Qnew, matrixSize);
    auto averageDuration = (double)duration.count();
    auto gFlops = flops / averageDuration;
    std::cout << "Matrix Multiplication GFLOPS: " << gFlops << std::endl;
    std::cout << "Matrix Multiplication FLOP: " << flops << std::endl;
    std::cout << "Matrix Multiplication Duration (s): " << averageDuration / 1e9 << std::endl;

    return 0;
}
