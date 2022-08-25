/* ====================================================
 * CUTEst interface for dune-nonlinopt
 * (Based on CUTEr/CUTEst interface files)
 * ====================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include<iomanip>

#include<dune/nonlinopt/nonlinopt.hh>
#include<dune/nonlinopt/vectorclass.hh>

#define MAXLINE 256

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "cutest.h"

  /* global variables */
  integer CUTEst_nvar;        /* number of variables */
  integer CUTEst_ncon;        /* number of constraints */

  class Problem
    : public Dune::NonlinOpt::ProblemBase<double,Dune::NonlinOpt::VectorClass<double>>
  {
    private:

      mutable int iter = 0;

    public:

      using Real = double;
      using Point = Dune::NonlinOpt::VectorClass<Real>;

      std::size_t dim() const override
      {
        return CUTEst_nvar;
      }

      Point zero() const override
      {
        return Point(CUTEst_nvar);
      }

      double value(const Point& x, bool subsequent = false) const override
      {
        double f;
        integer status;
        CUTEST_ufn(&status, &CUTEst_nvar, &(const_cast<Point&>(x)[0]), &f);
        if ((status == 1) || (status == 2))
        {
          printf("** CUTEst error, status = %d, aborting\n", status);
          exit(status);
        }

        return f;
      }

      void gradient(const Point& x, Point& g, bool subsequent = false) const override
      {
        integer status;
        if (g.size() != CUTEst_nvar)
          g = zero();
        CUTEST_ugr(&status, &CUTEst_nvar, &(const_cast<Point&>(x)[0]), &(g[0]));
        if ((status == 1) || (status == 2))
        {
          printf("** CUTEst error, status = %d, aborting\n", status);
          exit(status);
        }
      }

      double value_and_grad(const Point& x, Point& g) const override
      {
        double f;
        integer status;
        bool grad = true;
        if (g.size() != CUTEst_nvar)
          g = zero();
        CUTEST_uofg(&status, &CUTEst_nvar, &(const_cast<Point&>(x)[0]), &f, &(g[0]), &grad);

        if ((status == 1) || (status == 2))
        {
          printf("** CUTEst error, status = %d, aborting\n", status);
          exit(status);
        }

        return f;
      }

      void hook(std::size_t iter_, const Point& point,
          Real value, const Point& descent) const override
      {
        iter = iter_;
      }

      int iteration() const
      {
        return iter;
      }
  };

  using Real   = typename Problem::Real;
  using Point  = typename Problem::Point;
  using Solver = Dune::NonlinOpt::UnconstrainedOptimization<Real,Point>;

  int solve(Solver& solver, Problem& problem, double* x, int max_it)
  {
    int err = 0;

    Point x_vec(CUTEst_nvar);
    for (int i = 0; i < CUTEst_nvar; i++)
      x_vec[i] = x[i];

    std::cout << "starting iteration" << std::endl;
    try
    {
      solver.hard_reset(problem);
      int i = 0;
      while (i < max_it)
      {
        bool converged = solver.step(problem, x_vec);
        if (converged)
          break;
        i++;
      }
      if (i == max_it)
        err = 1;
      std::cout << "done at iteration " << i << std::endl;
    }
    catch(Dune::NonlinOpt::NoDescentDirection)
    {
      std::cout << "err: no descent" << std::endl;
      err = 251;
    }
    catch(Dune::NonlinOpt::BracketingFailed)
    {
      std::cout << "err: bracketing failed" << std::endl;
      err = 252;
    }
    catch(Dune::NonlinOpt::SlopeAlwaysNegative)
    {
      std::cout << "err: slope always negative" << std::endl;
      err = 253;
    }
    catch(Dune::NonlinOpt::LinesearchFailed)
    {
      std::cout << "err: linesearch failed" << std::endl;
      err = 254;
    }
    for (int i = 0; i < CUTEst_nvar; i++)
      x[i] = x_vec[i];

    std::cout << "done without errors" << std::endl;
    return err;
  }

  /* main program */
  int MAINENTRY(void)
  {
    const char* fname = "OUTSDIF.d"; /* CUTEst data file */
    integer funit = 42;              /* FORTRAN unit number for OUTSDIF.d */
    integer io_buffer = 11;          /* FORTRAN unit for internal i/o */
    integer iout = 6;                /* FORTRAN unit number for error output */
    integer ierr;                    /* Exit flag from OPEN and CLOSE */
    integer status;                  /* Exit flag from CUTEst tools */

    doublereal* x, * bl, * bu;
    char*       pname;
    logical     constrained = FALSE_;
    doublereal  calls[7], cpu[2];

    FILE* spec;

    /* Open problem description file OUTSDIF.d */
    ierr = 0;

    FORTRAN_open(&funit, fname, &ierr);
    if(ierr)
    {
      printf("Error opening file OUTSDIF.d.\nAborting.\n");
      exit(1);
    }

    /* Determine problem size */
    CUTEST_cdimen(&status, &funit, &CUTEst_nvar, &CUTEst_ncon);
    if (status)
    {
      printf("** CUTEst error, status = %d, aborting\n", status);
      exit(status);
    }
    /* Determine whether to call constrained or unconstrained tools */
    if(CUTEst_ncon)
      constrained = TRUE_;

    /* stop if the problem has constraints */
    if(constrained)
    {
      printf (" ** the problem %s has %i constraints\n",
          pname, CUTEst_ncon);
      printf ("    dune-nonlinopt is for unconstrained optimization\n");
      exit(-1);
    }

    /* Reserve memory for variables, bounds, and multipliers */
    /* and call appropriate initialization routine for CUTEst */
    MALLOC(x,  CUTEst_nvar, doublereal);
    MALLOC(bl, CUTEst_nvar, doublereal);
    MALLOC(bu, CUTEst_nvar, doublereal);
    CUTEST_usetup(&status, &funit, &iout, &io_buffer, &CUTEst_nvar,
        x, bl, bu);
    if (status)
    {
      printf("** CUTEst error, status = %d, aborting\n", status);
      exit(status);
    }

    /* Get problem name */
    MALLOC(pname, FSTRING_LEN+1, char);
    CUTEST_probname(&status, pname);
    if (status)
    {
      printf("** CUTEst error, status = %d, aborting\n", status);
      exit(status);
    }

    /* Make sure to null-terminate problem name */
    pname[FSTRING_LEN] = '\0';
    int i = FSTRING_LEN - 1;
    while(i-- > 0 && pname[i] == ' ')
    {
      pname[i] = '\0';
    }

    Solver solver(10000,10,1);

    Problem problem;
    /* Call the optimizer */
    int exit_code = solve(solver,problem,x,10000);

    /* Get CUTEst statistics */
    CUTEST_creport( &status, calls, cpu);
    if (status)
    {
      printf("** CUTEst error, status = %d, aborting\n", status);
      exit(status);
    }

    /* Evaluate final function value / gradient*/
    Point x_vec(CUTEst_nvar);
    for (int i = 0; i < CUTEst_nvar; i++)
      x_vec[i] = x[i];
    const double f_final = solver.get_value();
    const Point& g_final = solver.get_gradient();

    printf(" *********************** CUTEst statistics ************************\n");
    printf(" Code used               : dune-nonlinopt\n");
    printf(" Problem                 : %-s\n", pname);
    printf(" # variables             = %10d\n", CUTEst_nvar);
    printf(" # iterations            = %10d\n", problem.iteration());
    printf(" # objective functions   = %10.7g\n", calls[0]);
    printf(" # objective gradients   = %10.7g\n", calls[1]);
    printf(" # objective hessians    = %10.7g\n", calls[2]);
    printf(" # evaluations           = %10.7g\n", calls[0] + calls[1]);
    printf(" Exit code               = %10d\n", exit_code);
    printf(" Final f                 = %10.7E\n",f_final);
    printf(" Final ||g||             = %10.7E\n",std::sqrt(g_final*g_final));
    printf(" Set up time             = %10.7f seconds\n", cpu[0]);
    printf(" Solve time              = %10.7f seconds\n", cpu[1]);
    printf(" ******************************************************************\n");

    ierr = 0;
    FORTRAN_close(&funit, &ierr);
    if(ierr)
    {
      printf("Error closing %s on unit %d.\n", fname, funit);
      printf("Trying not to abort.\n");
    }

    /* Free workspace */
    FREE(pname);
    FREE(x); FREE(bl); FREE(bu);

    CUTEST_uterminate(&status);

    return 0;
  }

#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif
