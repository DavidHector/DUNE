// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ISTL_EIGENVALUE_TEST_MATRIXINFO_HH
#define DUNE_ISTL_EIGENVALUE_TEST_MATRIXINFO_HH

#include <cmath>    // provides std::abs and std::sqrt
#include <cassert>  // provides assert
#include <limits>
#include <iostream>  // provides std::cout, std::endl

#include <dune/common/exceptions.hh>  // provides DUNE_THROW(...), Dune::Exception
#include <dune/common/fvector.hh>     // provides Dune::FieldVector

#include <dune/istl/blocklevel.hh>       // provides Dune::blockLevel
#include <dune/istl/bvector.hh>          // provides Dune::BlockVector
#include <dune/istl/superlu.hh>          // provides Dune::SuperLU
#include <dune/istl/preconditioners.hh>  // provides Dune::SeqGS
#include <dune/istl/solvers.hh>          // provides Dune::BiCGSTABSolver
#include <dune/istl/matrixmatrix.hh>     // provides Dune::transposeMatMultMat(...)

#include "../poweriteration.hh"  // provides Dune::PowerIteration_Algorithms


/**
 * \brief Class template which yields information related to a square
 *        matrix like its spectral (i.e. 2-norm) condition number.
 *
 * \todo The current implementation is limited to DUNE-ISTL
 *       BCRSMatrix types with blocklevel 2. An extension to
 *       blocklevel >= 2 might be provided in a future version.
 *
 * \tparam BCRSMatrix Type of a DUNE-ISTL BCRSMatrix whose properties
 *                    shall be considered; is assumed to have blocklevel
 *                    2 with square blocks.
 *
 * \author Sebastian Westerheide.
 */
template <typename BCRSMatrix>
class MatrixInfo
{
public:
  //! Type of the underlying field of the matrix
  typedef typename BCRSMatrix::field_type Real;

public:
  /**
   * \brief Construct from required parameters.
   *
   * \param[in] m                       The DUNE-ISTL BCRSMatrix
   *                                    whose properties shall be
   *                                    considered; is assumed to
   *                                    be square.
   * \param[in] b                       The rigth hand side DUNE-ISTL BCRSMatrix
   * \param[in] verbose                 Verbosity setting.
   * \param[in] arppp_a_verbosity_level Verbosity setting of the
   *                                    underlying ARPACK++ algorithms.
   * \param[in] pia_verbosity_level     Verbosity setting of the
   *                                    underlying power iteration
   *                                    based algorithms.
   */
  MatrixInfo (const BCRSMatrix& m,
              const BCRSMatrix& b,
              const bool verbose = false,
              const unsigned int arppp_a_verbosity_level = 0,
              const unsigned int pia_verbosity_level = 0)
    : m_(m),
      b_(b),
      verbose_(verbose),
      arppp_a_verbosity_level_(arppp_a_verbosity_level*verbose),
      pia_verbosity_level_(pia_verbosity_level*verbose),
      cond_2_(-1.0), symmetricity_assumed_(false)
  {
    // assert that BCRSMatrix type has blocklevel 2
    static_assert
      (Dune::blockLevel<BCRSMatrix>() == 2,
       "Only BCRSMatrices with blocklevel 2 are supported.");

    // assert that BCRSMatrix type has square blocks
    static_assert
      (BCRSMatrix::block_type::rows == BCRSMatrix::block_type::cols,
       "Only BCRSMatrices with square blocks are supported.");

    // assert that m_ is square
    const int nrows = m_.M() * BCRSMatrix::block_type::rows;
    const int ncols = m_.N() * BCRSMatrix::block_type::cols;
    if (nrows != ncols)
      DUNE_THROW(Dune::Exception,"Matrix is not square ("
                 << nrows << "x" << ncols << ").");
  }

  //! Type of block vectors compatible with the rows of a BCRSMatrix
  //! object and its columns
  static const int bvBlockSize = BCRSMatrix::block_type::rows;
  typedef Dune::FieldVector<Real,bvBlockSize> BlockVectorBlock;
  typedef Dune::BlockVector<BlockVectorBlock> BlockVector;

  //! computes smallest EVs
  inline Real computeSmallestEV () const
  {
    // 1) allocate memory for largest and smallest magnitude eigenvalue
    //    as well as the spectral (i.e. 2-norm) condition number
    Real lambda_min_gen{}, lambda_min{}, cond_2{};

    // 2) allocate memory for starting vectors and approximated
    //    eigenvectors
    BlockVector x(m_.M());

    // 4) setup power iteration based iterative eigenvalue algorithms
    typedef Dune::PowerIteration_Algorithms<BCRSMatrix,BlockVector> PIA;
    PIA pia(m_,20000,pia_verbosity_level_);
    static const bool avoidLinSolverCrime = true;

#if HAVE_SUPERLU
    // 5) select a linear solver for power iteration based iterative
    //    eigenvalue algorithms
    typedef Dune::SuperLU<BCRSMatrix> PIALS;
    const unsigned int piaLS_verbosity = 0;
    PIALS piaLS(pia.getIterationMatrix(),piaLS_verbosity);
#else
    // 5) select a linear solver for power iteration based iterative
    //    eigenvalue algorithms
    typedef Dune::SeqGS<BCRSMatrix,
                        typename PIA::IterationOperator::domain_type,
                        typename PIA::IterationOperator::range_type> PIAPC;
    PIAPC piaPC(pia.getIterationMatrix(),2,1.0);
    const double piaLS_reduction = 1e-02;
    const unsigned int piaLS_max_iter = 1000;
    const unsigned int piaLS_verbosity = 0;
    typedef Dune::BiCGSTABSolver<typename PIA::IterationOperator::domain_type> PIALS;
    PIALS piaLS(pia.getIterationOperator(),piaPC,
                piaLS_reduction,piaLS_max_iter,piaLS_verbosity);
#endif  // HAVE_SUPERLU

    const Real epsilonPI    = 1e-02;
    const Real gamma = 0.0;
    x = 1.0;
    // 6.1) perform power iteration for smallest magnitude
    //       eigenvalue (predictor)
    pia.applyInverseIteration(gamma, epsilonPI,piaLS, x,lambda_min);

    // 10) output smallest magnitude eigenvalue
    if (verbose_)
      std::cout << "    Smallest magnitude eigenvalue (standard implementation) λ_min = "
                << lambda_min << std::endl;

    x = 1.0;
    // 6.1) perform power iteration for smallest generalised magnitude
    //       eigenvalue (predictor)
    pia.applyGeneralInverseIteration(b_, gamma, epsilonPI,piaLS, x,lambda_min_gen);

    // 10) output smallest magnitude eigenvalue
    if (verbose_)
      std::cout << "    Smallest magnitude eigenvalue (my implementation) λ_min_gen = "
                << lambda_min_gen << std::endl;
    return lambda_min;
  }


protected:
  // parameters related to computation of matrix information
  const BCRSMatrix& m_;
  const BCRSMatrix& b_;

  // verbosity setting
  const bool verbose_;
  const unsigned int arppp_a_verbosity_level_;
  const unsigned int pia_verbosity_level_;

  // memory for storing matrix information
  // (mutable as matrix information is computed on demand)
  mutable Real cond_2_;
  mutable bool symmetricity_assumed_;
};


#endif  // DUNE_ISTL_EIGENVALUE_TEST_MATRIXINFO_HH
