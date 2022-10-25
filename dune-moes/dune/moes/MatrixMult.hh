#ifndef DUNE_MOES_MATRIXMULT_HH
#define DUNE_MOES_MATRIXMULT_HH

#include <dune/moes/vectorclass/vectorclass.h>
#include <dune/istl/bcrsmatrix.hh>
#include <iostream>

/*
    input: M: Matrix to calculate EVs for, Q: 

*/
template <typename MT>
void MultQ(MT &M, Vec4d *Qin, Vec4d *Qout, size_t qCols, size_t N)
{
    // How to iterate over M?
    // So I think for a BCRS Matrix, this only iterates over the blocks (which are fieldmatrices I guess)
    // Set Qout to zero
    for (size_t qCol = 0; qCol < qCols; qCol++)
    {
        for (size_t qRow = 0; qRow < 2 * N; qRow++)
        {
            Qout[qCol * 2 * N + qRow] = 0.0;
        }
    }

    int rows = 0;
    auto endRow = M.end();
    for (auto rowIterator = M.begin(); rowIterator != endRow; rowIterator++)
    {
        auto endCol = (*rowIterator).end();
        auto mBrI = rowIterator.index(); // Matrix Block Row Index
        rows++;
        int cols = 0;
        for (auto colIterator = (*rowIterator).begin(); colIterator != endCol; colIterator++)
        {
            cols++;
            // *colIterator is probably the FieldMatrix (which uses the dense matrix multAssign)
            auto mBcI = rowIterator.index();
            auto fMatRows = (*colIterator).rows; // Should be blockSize of the underlying field matrix
            auto fMatCols = (*colIterator).cols;
            for (auto i = 0; i < fMatRows; i++)
            {
                auto mIr = mBrI * fMatRows + i;
                for (size_t qCol = 0; qCol < qCols; qCol++)
                {
                    for (auto j = 0; j < fMatCols; j++)
                    {
                        auto mIc = mBcI * fMatCols + j;
                        auto qinIndex = qCol * 2 * N + mIc * 2;  // columns in the matrix are rows in the vector
                        auto qoutIndex = qCol * 2 * N + mIr * 2; //col row
                        Qout[qoutIndex] += (*colIterator)[i][j] * Qin[qinIndex];
                        Qout[qoutIndex + 1] += (*colIterator)[i][j] * Qin[qinIndex + 1];
                        // Problem: Q has BlockSize 2, so at every qrow (i.e. mIr/mIc, I have to also look at the other block)
                        // basically what I want is:
                        // Qout[qCol][mIr] += M[mBrI][mBcI][i][j] * Qin[qCol][mIc];
                        // reorder the loop to avoid striding

                        // (*colIterator)[i][j] should be the actual entry, while rowIterator.index() amd colIterator.index() are the blockIndices of the bcrsmatrix
                        // How do I do the multiplication?
                        // We are going going through the columns, but are always fMatRows wide. I have to multiply it with everything anyway, so maybe the rows dont matter
                        // What to do about the fact that the matrix has doubles as entries. But Q has Vec4ds
                    }
                }
            }
        }
    }
}
/*
    Matrix multiplication with double Q 
*/

template <typename MT>
void MultQ(MT &M, const double *Qin, double *Qout, size_t qCols, size_t N)
{
    Vec4d v, u;
    Vec4d zeroVec = 0.0;
    Vec4d productFirst, productSecond;
    Vec4d inFirst, inSecond;
    Vec4d outFirst, outSecond;
    Vec4d entryM;
    // How to iterate over M?
    // So I think for a BCRS Matrix, this only iterates over the blocks (which are fieldmatrices I guess)
    // Set Qout to zero
    for (size_t i = 0; i < qCols * N * 8; i += 4)
    {
        zeroVec.store(&Qout[i]);
    }

    int rows = 0;
    auto endRow = M.end();
    for (auto rowIterator = M.begin(); rowIterator != endRow; rowIterator++)
    {
        auto endCol = (*rowIterator).end();
        auto mBrI = rowIterator.index(); // Matrix Block Row Index
        rows++;
        int cols = 0;
        for (auto colIterator = (*rowIterator).begin(); colIterator != endCol; colIterator++)
        {
            cols++;
            // *colIterator is probably the FieldMatrix (which uses the dense matrix multAssign)
            auto mBcI = rowIterator.index();     // Matrix Block Column Index
            auto fMatRows = (*colIterator).rows; // Should be blockSize of the underlying field matrix
            auto fMatCols = (*colIterator).cols;
            // std::cout << "fMatRows = " << fMatRows << ", fMatCols = " << fMatCols << std::endl;
            for (auto i = 0; i < fMatRows; i++)
            {
                auto mIr = mBrI * fMatRows + i; // Matrix Row Index
                for (size_t qCol = 0; qCol < qCols; qCol++)
                {
                    productFirst = 0.0;
                    productSecond = 0.0;
                    auto qoutIndex = qCol * N * 8 + mIr * 8; //col row  this should go to 64, but only goes to 24, why?
                    //std::cout << "qoutIndex = " << qoutIndex << ", mIr = " << mIr << ", mBrI = " << mBrI << std::endl;
                    for (auto j = 0; j < fMatCols; j++)
                    {
                        auto mIc = mBcI * fMatCols + j; // Matrix Column Index, this is somehow very wrong
                        // std::cout << "qCol = " << qCol << ", mBcI = " << mBcI << ", fMatCols = " << fMatCols << ", mIc = " << mIc << ", mIr = " << mIr << std::endl;
                        auto qinIndex = qCol * N * 8 + mIc * 8; // columns in the matrix are rows in the vector
                        inFirst.load(&Qin[qinIndex]);
                        inSecond.load(&Qin[qinIndex + 4]);
                        entryM = (*colIterator)[i][j]; // This seems fine
                        //std::cout << entryM << std::endl;

                        productFirst += entryM * inFirst;
                        productSecond += entryM * inSecond;
                    }
                    productFirst.store(&Qout[qoutIndex]);
                    productSecond.store(&Qout[qoutIndex + 4]);
                }
            }
        }
    }
}

#endif