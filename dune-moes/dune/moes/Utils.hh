#ifndef DUNE_MOES_UTILS_HH
#define DUNE_MOES_UTILS_HH

#include <dune/moes/vectorclass/vectorclass.h>
/**
 * @brief Solves L x = b
 * 
 * @param b Multivector to go in (Qin)
 * @param x Multivector to go out (Qout)
 * @param Lp 
 * @param Lj 
 * @param Lx 
 */
void solveL(const std::shared_ptr<double[]> &b, std::shared_ptr<double[]> &x, const int64_t *Lp, const int64_t *Lj, const double *Lx, const int64_t n_row, const size_t qCols)
{
    // Dont have to initialize x, because every value will be written
    const size_t qSize = n_row * qCols * 8;
    Vec4d LEntry, bFirst, bSecond, xFirst, xSecond;
    size_t bIndex, xIndex, LCol;
    // iterate over the rows of L
    for (int64_t row = 0; row < n_row; row++)
    {
        // Not really sure what is in Lp -> Index ranges for the values in the row
        // std::cout << "Lp[" << row << "] = " << Lp[row] << std::endl;
        // Multiply matrix rows with columns of Q
        for (size_t qCol = 0; qCol < qCols; qCol++)
        {
            xFirst = 0.0;
            xSecond = 0.0;
            bIndex = 8 * qCol * n_row + 8 * row;
            // Go through the previous (solved entries) until col = row, where the right side comes into play
            // If we can trust (have to check) that columns are in order, then we can do the back insertion via right-side writing
            for (int64_t j = Lp[row]; j < Lp[row + 1] - 1; j++) // have to check if index starts with 1 or 0
            {
                LCol = Lj[j];
                LEntry = Lx[j]; // broadcast
                xIndex = 8 * qCol * n_row + 8 * LCol;
                bFirst.load(&x[xIndex]);
                bSecond.load(&x[xIndex + 4]);
                xFirst -= LEntry * bFirst; // This here is for the non-diagonal elements that have to use the previous solutions
                xSecond -= LEntry * bSecond;
            }
            // Diagonal entry
            LCol = Lj[Lp[row + 1] - 1];   // this is for the entry corresponding right side b (diagonal element)
            LEntry = Lx[Lp[row + 1] - 1]; // This should always be 1 (Since L is a unit lower triangular matrix)
            // Checked this: It does work as intended, no non-1 diagonals were found that could cause the nans
            /*
            if (LEntry[0] != 1.0)
            {
                std::cout << "solveL: Encountered non-1 diagonal: " << LEntry[0] << std::endl;
            }
            */
            bFirst.load(&b[bIndex]);
            bSecond.load(&b[bIndex + 4]);
            xFirst += bFirst;
            xFirst /= LEntry;
            xSecond += bSecond;
            xSecond /= LEntry;
            xFirst.store(&x[bIndex]);
            xSecond.store(&x[bIndex + 4]);
        }
    }
}

/**
 * @brief Solves U x = b, should basically be like solveL, just in reverse, U is in CCS (Compressed Column Storage) format.
 * 
 * @param b Multivector to go in 
 * @param x Multivector to go out
 * @param Up 
 * @param Ui 
 * @param Ux 
 */
void solveU(const std::shared_ptr<double[]> &b, std::shared_ptr<double[]> &x, const int64_t *Up, const int64_t *Ui, const double *Ux, const int64_t n_col, const size_t qCols)
{
    // initialize x
    for (size_t i = 0; i < 8 * n_col * qCols; i++)
    {
        x[i] = 0.0;
    }
    Vec4d xFirst, xSecond, UEntry, bFirst, bSecond, tmpFirst, tmpSecond;
    size_t bIndex, xIndex, URow;
    for (size_t qCol = 0; qCol < qCols; qCol++)
    {
        // Go from right to left and bottom to top
        for (int64_t col = n_col - 1; col >= 0; col--)
        {
            // Go from bottom to top
            // Diagonal Entry
            URow = Ui[Up[col + 1] - 1];
            UEntry = Ux[Up[col + 1] - 1];
            // Checked this: It does work as intended, no small diagonals were found that could cause the nans
            /*
            if (std::abs(UEntry[0]) < 0.00005)
            {
                std::cout << "solveU: Encountered very small diagonal: " << UEntry[0] << std::endl;
            }
            */
            bIndex = 8 * n_col * qCol + 8 * URow;
            bFirst.load(&b[bIndex]);
            bSecond.load(&b[bIndex + 4]);
            xFirst.load(&x[bIndex]);
            xSecond.load(&x[bIndex + 4]);
            xFirst += bFirst;
            xSecond += bSecond;
            xFirst /= UEntry;
            xSecond /= UEntry;
            xFirst.store(&x[bIndex]);
            xSecond.store(&x[bIndex + 4]);
            for (int64_t i = Up[col + 1] - 2; i >= Up[col]; i--)
            {
                URow = Ui[i];
                UEntry = Ux[i];
                xIndex = 8 * n_col * qCol + 8 * URow;
                tmpFirst.load(&x[xIndex]);
                tmpSecond.load(&x[xIndex + 4]);
                tmpFirst -= UEntry * xFirst;
                tmpSecond -= UEntry * xSecond;
                tmpFirst.store(&x[xIndex]);
                tmpSecond.store(&x[xIndex + 4]);
            }
        }
    }
}

/**
 * @brief Applies the pivot x = Pb with P being a vector of row indices as provided by Umfpack
 * 
 * @param b Mutlivector to go in
 * @param x Multivector to go out
 * @param P Pivot vector
 * @param N Length of the Pivot vector/Multivectors
 * @param qCols Width of the Multivectors
 */
void applyPivotP(const std::shared_ptr<double[]> &b, std::shared_ptr<double[]> &x, const int64_t *P, const size_t N, const size_t qCols)
{
    Vec4d bFirst, bSecond;
    size_t bIndex, xIndex;
    for (size_t qCol = 0; qCol < qCols; qCol++)
    {
        for (size_t i = 0; i < N; i++)
        {
            /* load old vector at correct indices and write into new vector at i */
            bIndex = 8 * qCol * N + 8 * P[i];
            xIndex = 8 * qCol * N + 8 * i; // This can be much optimized
            bFirst.load(&b[bIndex]);
            bSecond.load(&b[bIndex + 4]);
            bFirst.store(&x[xIndex]);
            bSecond.store(&x[xIndex + 4]);
        }
    }
}

/**
 * @brief Applies the pivot x = Qb with Q being a vector of column indices as provided by Umfpack
 * 
 * @param b Mutlivector to go in
 * @param x Multivector to go out
 * @param Q Pivot vector
 * @param N Length of the Pivot vector/Multivectors
 * @param qCols Width of the Multivectors
 */

void applyPivotQ(const std::shared_ptr<double[]> &b, std::shared_ptr<double[]> &x, const int64_t *Q, const size_t N, const size_t qCols)
{
    Vec4d bFirst, bSecond;
    size_t bIndex, xIndex;
    for (size_t qCol = 0; qCol < qCols; qCol++)
    {
        for (size_t i = 0; i < N; i++)
        {
            xIndex = 8 * qCol * N + 8 * Q[i];
            bIndex = 8 * qCol * N + 8 * i; // This can be much optimized
            bFirst.load(&b[bIndex]);
            bSecond.load(&b[bIndex + 4]);
            bFirst.store(&x[xIndex]);
            bSecond.store(&x[xIndex + 4]);
        }
    }
}

/**
 * @brief Normalizes x
 * 
 * @param x Multivector to normalize
 * @param N Length of the Multivector
 * @param qCols rhsWidth/8 (number of columns)
 */
void normalize(std::shared_ptr<double[]> &x, const size_t N, const size_t qCols)
{
    Vec4d xFirst, xSecond, normFirst, normSecond, oneVec;
    size_t xIndex;
    oneVec = 1.0; // broadcast for efficent division later
    for (size_t qCol = 0; qCol < qCols; qCol++)
    {
        normFirst = 0.0;
        normSecond = 0.0;
        for (size_t i = 0; i < N; i++)
        {
            xIndex = 8 * qCol * N + 8 * i;
            xFirst.load(&x[xIndex]);
            xSecond.load(&x[xIndex + 4]);
            normFirst += square(xFirst);
            normSecond += square(xSecond);
        }
        normFirst = oneVec / sqrt(normFirst);
        normSecond = oneVec / sqrt(normSecond);
        for (size_t i = 0; i < N; i++)
        {
            xIndex = 8 * qCol * N + 8 * i;
            xFirst.load(&x[xIndex]);
            xSecond.load(&x[xIndex + 4]);
            xFirst *= normFirst;
            xSecond *= normSecond;
            xFirst.store(&x[xIndex]);
            xSecond.store(&x[xIndex + 4]);
        }
    }
}
void applyScalingMult(const std::shared_ptr<double[]> &b, std::shared_ptr<double[]> &x, const double *Rs, const size_t N, const size_t qCols)
{
    Vec4d REntry, xFirst, xSecond;
    size_t xIndex = 0;
    for (size_t qCol = 0; qCol < qCols; qCol++)
    {
        for (size_t row = 0; row < N; row++)
        {
            REntry = Rs[row]; // broadcast
            // xIndex = 8 * qCol * N + 8 * row;
            xFirst.load(&b[xIndex]);
            xSecond.load(&b[xIndex + 4]);
            xFirst *= REntry;
            xSecond *= REntry;
            xFirst.store(&x[xIndex]);
            xSecond.store(&x[xIndex + 4]);
            xIndex += 8;
        }
    }
};
void applyScalingInverse(const std::shared_ptr<double[]> &b, std::shared_ptr<double[]> &x, const double *Rs, const size_t N, const size_t qCols)
{
    Vec4d REntry, xFirst, xSecond;
    size_t xIndex = 0;
    for (size_t qCol = 0; qCol < qCols; qCol++)
    {
        for (size_t row = 0; row < N; row++)
        {
            REntry = Rs[row]; // broadcast
            // xIndex = 8 * qCol * N + 8 * row;
            xFirst.load(&b[xIndex]);
            xSecond.load(&b[xIndex + 4]);
            xFirst /= REntry;
            xSecond /= REntry;
            xFirst.store(&x[xIndex]);
            xSecond.store(&x[xIndex + 4]);
            xIndex += 8;
        }
    }
};

/**
 * @brief Applies the scaling vector R, which represents a diagonal matrix
 * 
 * @param b 
 * @param x 
 * @param Rs 
 * @param N 
 * @param qCols 
 * @param do_recip 
 */
void applyScaling(const std::shared_ptr<double[]> &b, std::shared_ptr<double[]> &x, const double *Rs, const size_t N, const size_t qCols, const int do_recip)
{
    if (do_recip == 1)
    {
        applyScalingMult(b, x, Rs, N, qCols);
    }
    else
    {
        applyScalingInverse(b, x, Rs, N, qCols);
    }
}

#endif