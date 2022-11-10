#ifndef DUNE_MOES_HH
#define DUNE_MOES_HH

#include <dune/moes/MatrixMult.hh>
#include <dune/moes/qrcol.hh>
#include <dune/moes/vectorclass/vectorclass.h>
#include <vector>

/*
    Checks, whether the sum of absolute differences between vector elements between iterations is lower than the given tolerance.
    
    Returns true if the check has passed, else returns false
*/
bool checkIterationTolerance(double *Qin, double *Qout, const size_t N, const size_t qCols, const double tolerance)
{
    double difference[8];
    Vec4d dOne, dTwo;
    Vec4d QinOne, QinTwo, QoutOne, QoutTwo;
    size_t qIndex;
    for (size_t col = 0; col < qCols; col++)
    {
        // reset difference
        dOne = 0.0;
        dTwo = 0.0;
        for (size_t row = 0; row < N; row++)
        {
            qIndex = col * N * 8 + row * 8;

            // loading both blocks
            QinOne.load(&Qin[qIndex]);
            QinTwo.load(&Qin[qIndex + 4]);
            QoutOne.load(&Qout[qIndex]);
            QoutTwo.load(&Qout[qIndex + 4]);

            dOne += QoutOne - QinOne;
            dTwo += QoutTwo - QinTwo;
        }
        dOne.store(&difference[0]);
        dTwo.store(&difference[4]);
        for (size_t i = 0; i < 4; i++)
        {
            if (std::abs(dOne[i]) > tolerance)
            {
                return false;
            }

            if (std::abs(dTwo[i]) > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

template <typename T>
void pointerSwap(T **a, T **b)
{
    T *tmp = *a;
    *a = *b;
    *b = tmp;
}

// add your classes here
template <typename MT>
void largestEVs(const MT &M, double *Q, const size_t qCols, const size_t N, const double tolerance, const size_t qrFrequency)
{
    bool stop = false;
    size_t iterationCounter = 0;
    size_t matrixSize = N * qCols * 8; // watch out, overflow error might occur, 8 because of col width
    double *Qtmp = new double[matrixSize];
    fillMatrixRandom(Qtmp, matrixSize);
    // printMatrix(Qtmp, N, qCols * 8);
    MultQ(M, Qtmp, Q, qCols, N);
    // printMatrix(Q, N, qCols * 8);
    while (!stop)
    {
        // TODO:
        // 1. Matrix multiplication
        // 2. call qr when iterationCounter%qrFrequency = 0;
        // 3. Check for tolerance between iterations (maybe some difference metric, after normalization (or qr call)).
        // 4. Switch Qin and Qout;

        // Matrix Multiplication
        // std::cout << "Before: Q[0] = " << Q[0] << std::endl;
        //std::cout << "Before: Qtmp[0] = " << Qtmp[0] << std::endl;
        MultQ(M, Qtmp, Q, qCols, N);
        // Why do the pointer swap, I could just do two multiplications in each step
        MultQ(M, Q, Qtmp, qCols, N);
        //std::cout << "After: Q[0] = " << Q[0] << std::endl;
        //std::cout << "After: Qtmp[0] = " << Qtmp[0] << std::endl;
        // Problem: The value of the entries is rising until it reaches nan, without converging -> there is probably something wrong with the multiplication
        // The matrix is only have filled by the end (WHY?)

        // Call QR Algorithm and check tolerance
        if (iterationCounter % qrFrequency == 0)
        {
            //std::cout << "Before QR: Q[0] = " << Q[0] << std::endl;
            //printMatrix(Q, N, qCols * 8);
            qrFixedBlockOptimizedDouble(Q, N, qCols, 2, 1);
            //std::cout << "After QR: Q[0] = " << Q[0] << std::endl;
            printMatrix(Q, N, qCols * 8);
            stop = checkIterationTolerance(Q, Qtmp, N, qCols, tolerance);
            if (stop)
            {
                std::cout << "largestEVs: Returning Q = " << std::endl;
                printMatrix(Q, N, qCols * 8);
                delete Qtmp;
                std::cout << "largestEVs took " << iterationCounter << " iterations to complete" << std::endl;
                return;
            }
        }
        iterationCounter++;
    }
    delete[] Qtmp;
}

template <typename MT>
void largestEVsIterative(const MT &M, double *Q, const size_t qCols, const size_t N, const size_t iterations, const size_t qrFrequency)
{
    bool stop = false;
    size_t matrixSize = N * qCols * 8; // watch out, overflow error might occur, 8 because of col width
    double *Qtmp = new double[matrixSize];
    fillMatrixRandom(Qtmp, matrixSize);
    // printMatrix(Qtmp, N, qCols * 8);
    MultQ(M, Qtmp, Q, qCols, N);
    // printMatrix(Q, N, qCols * 8);
    for (size_t i = 0; i < iterations; i++)
    {
        // TODO:
        // 1. Matrix multiplication
        // 2. call qr when iterationCounter%qrFrequency = 0;
        // 3. Check for tolerance between iterations (maybe some difference metric, after normalization (or qr call)).
        // 4. Switch Qin and Qout;

        // Matrix Multiplication
        // std::cout << "Before: Q[0] = " << Q[0] << std::endl;
        //std::cout << "Before: Qtmp[0] = " << Qtmp[0] << std::endl;
        MultQ(M, Qtmp, Q, qCols, N);
        // Why do the pointer swap, I could just do two multiplications in each step
        MultQ(M, Q, Qtmp, qCols, N);
        //std::cout << "After: Q[0] = " << Q[0] << std::endl;
        //std::cout << "After: Qtmp[0] = " << Qtmp[0] << std::endl;
        // Problem: The value of the entries is rising until it reaches nan, without converging -> there is probably something wrong with the multiplication
        // The matrix is only have filled by the end (WHY?)

        // Call QR Algorithm and check tolerance
        if (i % qrFrequency == 0)
        {
            //std::cout << "Before QR: Q[0] = " << Q[0] << std::endl;
            //printMatrix(Q, N, qCols * 8);
            qrFixedBlockOptimizedDouble(Q, N, qCols, 2, 1);
            //std::cout << "After QR: Q[0] = " << Q[0] << std::endl;
            printMatrix(Q, N, qCols * 8);
        }
    }
    std::cout << "largestEVs: Returning Q = " << std::endl;
    printMatrix(Q, N, qCols * 8);
    delete[] Qtmp;
    return;
}

template <typename MT>
void getEigenvalues(const MT &M, double *Q, const size_t qCols, const size_t N, std::vector<double> &EVs)
{
    size_t matrixSize = N * qCols * 8; // watch out, overflow error might occur, 8 because of col width
    double *Qtmp = new double[matrixSize];
    double EV = 0.0;
    size_t col, uIndex, offset;
    MultQ(M, Q, Qtmp, qCols, N);
    for (size_t EVIndex = 0; EVIndex < EVs.size(); EVIndex++)
    {
        col = EVIndex / 8;
        uIndex = col * N * 8 + offset;
        for (size_t n = 0; n < N; n++)
        {
            EV += Qtmp[uIndex] / Q[uIndex];
            uIndex += 8;
        }
        EVs[EVIndex] = EV / N;
    }
    delete[] Qtmp;
}
#endif // DUNE_MOES_HH
