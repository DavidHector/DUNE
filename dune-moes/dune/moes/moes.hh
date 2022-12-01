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
bool checkIterationTolerance(std::unique_ptr<double[]> &Qin, std::unique_ptr<double[]> &Qout, const size_t N, const size_t qCols, const double tolerance)
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
void largestEVs(const MT &M, std::unique_ptr<double[]> &Q, const size_t qCols, const size_t N, const double tolerance, const size_t qrFrequency)
{
    bool stop = false;
    size_t iterationCounter = 0;
    size_t matrixSize = N * qCols * 8; // watch out, overflow error might occur, 8 because of col width
    // double *Qtmp = new double[matrixSize];
    std::unique_ptr<double[]> Qtmp(new double[matrixSize]);
    fillMatrixRandom(Qtmp, matrixSize);
    // printMatrix(Qtmp, N, qCols * 8);
    fillMatrixRandom(M, Qtmp, Q, qCols, N);
    // printMatrix(Q, N, qCols * 8);
    while (!stop)
    {
        fillMatrixRandom(M, Qtmp, Q, qCols, N);
        // Why do the pointer swap, I could just do two multiplications in each step
        fillMatrixRandom(M, Q, Qtmp, qCols, N);

        // Call QR Algorithm and check tolerance
        if (iterationCounter % qrFrequency == 0)
        {
            //std::cout << "Before QR: Q[0] = " << Q[0] << std::endl;
            //printMatrix(Q, N, qCols * 8);
            qrFixedBlockOptimizedDouble(Q, N, qCols, 2, 1);
            //std::cout << "After QR: Q[0] = " << Q[0] << std::endl;
            // printMatrix(Q, N, qCols * 8);
            stop = checkIterationTolerance(Q, Qtmp, N, qCols, tolerance);
            if (stop)
            {
                std::cout << "largestEVs: Returning Q = " << std::endl;
                // printMatrix(Q, N, qCols * 8);
                // delete Qtmp;
                std::cout << "largestEVs took " << iterationCounter << " iterations to complete" << std::endl;
                return;
            }
        }
        iterationCounter++;
    }
}

template <typename MT>
void largestEVsIterative(const MT &M, std::unique_ptr<double[]> &Q, const size_t qCols, const size_t N, const size_t iterations, const size_t qrFrequency)
{
    bool stop = false;
    size_t matrixSize = N * qCols * 8; // watch out, overflow error might occur, 8 because of col width
    std::unique_ptr<double[]> Qtmp(new double[matrixSize]);
    fillMatrixRandom(Qtmp, matrixSize);
    powerIteration(M, Qtmp, Q, qCols, N);
    for (size_t i = 0; i < iterations; i++)
    {
        powerIteration(M, Qtmp, Q, qCols, N);
        // Why do the pointer swap, I could just do two multiplications in each step
        powerIteration(M, Q, Qtmp, qCols, N);

        // Call QR Algorithm and check tolerance
        if (i % qrFrequency == 0)
        {
            qrFixedBlockOptimizedDouble(Q, N, qCols, 2, 1);
        }
    }
    return;
}

template <typename MT>
void smallestEVsIterative(const MT &M, std::unique_ptr<double[]> &Q, const size_t qCols, const size_t N, const size_t iterations, const size_t qrFrequency)
{
    bool stop = false;
    size_t matrixSize = N * qCols * 8; // watch out, overflow error might occur, 8 because of col width
    std::unique_ptr<double[]> Qtmp(new double[matrixSize]);
    fillMatrixRandom(Qtmp, matrixSize);
    inversePowerIteration(M, Qtmp, Q, qCols, N);
    for (size_t i = 0; i < iterations; i++)
    {
        inversePowerIteration(M, Qtmp, Q, qCols, N);
        // Why do the pointer swap, I could just do two multiplications in each step
        inversePowerIteration(M, Q, Qtmp, qCols, N);

        // Call QR Algorithm and check tolerance
        if (i % qrFrequency == 0)
        {
            qrFixedBlockOptimizedDouble(Q, N, qCols, 2, 1);
        }
    }
    return;
}

/*
    getEigenvalues
    gets EVs by using the \mu_k = \frac{b_k^\ast A b_k}{b_k^\ast b_k} Rayleigh-Quotient
*/
template <typename MT>
void getEigenvalues(const MT &M, std::unique_ptr<double[]> &Q, const size_t qCols, const size_t N, std::vector<double> &EVs)
{
    size_t matrixSize = N * qCols * 8; // watch out, overflow error might occur, 8 because of col width
    std::unique_ptr<double[]> Qtmp(new double[matrixSize]);
    double EV = 0.0;
    double bkbk, u; // b_k^\ast b_k
    size_t col, uIndex, offset;
    MultQSimple(M, Q, Qtmp, qCols, N); // Qtmp = A b_k
    for (size_t EVIndex = 0; EVIndex < EVs.size(); EVIndex++)
    {
        col = EVIndex / 8;
        offset = EVIndex % 8;
        uIndex = col * N * 8 + offset;
        EV = 0.0;
        bkbk = 0.0;
        for (size_t qRow = 0; qRow < N; qRow++)
        {
            u = Q[uIndex];
            EV += u * Qtmp[uIndex];
            bkbk += u * u;
            uIndex += 8;
        }
        EVs[EVIndex] = EV / bkbk;
    }
}
#endif // DUNE_MOES_HH
