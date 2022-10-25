#ifndef DUNE_MOES_HH
#define DUNE_MOES_HH

#include <dune/moes/MatrixMult.hh>
#include <dune/moes/qrcol.hh>
#include <dune/moes/vectorclass/vectorclass.h>

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
    bool continue = true;
    size_t iterationCounter = 0;
    size_t matrixSize = N * qCols * 8; // watch out, overflow error might occur, 8 because of col width
    double *Qin = new double[matrixSize];
    double *Qout = Q;
    double *tmp;
    fillMatrixRandom(Qin, matrixSize);
    while (continue)
    {
        // TODO:
        // 1. Matrix multiplication
        // 2. call qr when iterationCounter%qrFrequency = 0;
        // 3. Check for tolerance between iterations (maybe some difference metric, after normalization (or qr call)).
        // 4. Switch Qin and Qout;

        // Matrix Multiplication
        MultQ(M, Qin, Qout, qCols, N);

        // Call QR Algorithm and check tolerance
        if (iterationCounter % qrFrequency == 0)
        {
            qrFixedBlockOptimizedDouble(Qout, N, qCols, 2, 1);
            continue = !checkIterationTolerance(Qin, Qout, N, qCols, tolerance);
            if (!continue && (Q == Qout))
            {
                return;
            }
        }
        // Switcheroo
        pointerSwap(&Qin, &Qout);
        iterationCounter++;
    }
}
#endif // DUNE_MOES_HH
