#include <config.h>

#include <iostream>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/istl/bvector.hh>

#include "../../test/laplacian.hh"
#include "../../test/identity.hh"
#include "matrixinfo_general.hh"


int main (int argc, char** argv)
{
  try
  {
    typedef double FIELD_TYPE;

    static const int BS = 1;
    std::size_t N = 60;

    if (argc > 1)
      N = atoi(argv[1]);
    std::cout << "testing for N = " << N << ", BS = " << BS << std::endl;

    typedef Dune::FieldMatrix<FIELD_TYPE,BS,BS> MatrixBlock;
    typedef Dune::BCRSMatrix<MatrixBlock> BCRSMat;

    BCRSMat mat;
    setupLaplacian(mat,N);
    BCRSMat id;
    setupIdentity(id, N);

    Dune::Timer watch;

    const bool verbose = true;
    const unsigned int arppp_a_verbosity_level = 2;
    const unsigned int pia_verbosity_level = 1;
    MatrixInfo<BCRSMat> matrixInfo
      (mat, id, verbose,arppp_a_verbosity_level,pia_verbosity_level);

    watch.reset();
    matrixInfo.computeSmallestEV();
    std::cout << "computation of condition number took " << watch.elapsed()
              <<" seconds" << std::endl;

    return 0;
  }
  catch (std::exception& e)
  {
    std::cout << "ERROR: " << e.what() << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
