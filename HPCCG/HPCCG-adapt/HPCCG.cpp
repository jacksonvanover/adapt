
//@HEADER
// ************************************************************************
// 
//               HPCCG: Simple Conjugate Gradient Benchmark Code
//                 Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/////////////////////////////////////////////////////////////////////////

// Routine to compute an approximate solution to Ax = b where:

// A - known matrix stored as an HPC_Sparse_Matrix struct

// b - known right hand side vector

// x - On entry is initial guess, on exit new approximate solution

// max_iter - Maximum number of iterations to perform, even if
//            tolerance is not met.

// tolerance - Stop and assert convergence if norm of residual is <=
//             to tolerance.

// niters - On output, the number of iterations actually performed.

/////////////////////////////////////////////////////////////////////////

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <cmath>
#include "mytimer.hpp"
#include "HPCCG.hpp"
#include "adapt-impl.cpp"

#define TICK()  t0 = mytimer() // Use TICK and TOCK to time a code section
#define TOCK(t) t += mytimer() - t0

int HPCCG(HPC_Sparse_Matrix * A,
	  const AD_real * const b, AD_real * const x,
	  const int max_iter, const double tolerance, int &niters, AD_real & normr,
	  double * times)

{
  double t_begin = mytimer();  // Start timing right away

  double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;
#ifdef USING_MPI
  double t5 = 0.0;
#endif
  int nrow = A->local_nrow;
  int ncol = A->local_ncol;

  AD_real * r = new AD_real [nrow];
  AD_real * p = new AD_real [ncol]; // In parallel case, A is rectangular
  AD_real * Ap = new AD_real [nrow];

  normr = 0.0;
  AD_real rtrans = 0.0;
  AD_real oldrtrans = 0.0;

#ifdef USING_MPI
  int rank; // Number of MPI processes, My process ID
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  int rank = 0; // Serial case (not using MPI)
#endif

  int print_freq = max_iter/10; 
  if (print_freq>50) print_freq=50;
  if (print_freq<1)  print_freq=1;

  // p is of length ncols, copy x to p for sparse MV operation
  TICK(); waxpby(nrow, 1.0, x, 0.0, x, p); TOCK(t2);
#ifdef USING_MPI
  TICK(); exchange_externals(A,p); TOCK(t5); 
#endif
  TICK(); HPC_sparsemv(A, p, Ap); TOCK(t3);
  TICK(); waxpby(nrow, 1.0, b, -1.0, Ap, r); TOCK(t2);
  TICK(); ddot(nrow, r, r, &rtrans, t4); TOCK(t1);
  normr = sqrt(AD_value(rtrans));

  if (rank==0) cout << "Initial Residual = "<< AD_value(normr) << endl;

  for(int k=1; k<max_iter && normr > tolerance; k++ )
    {
      if (k == 1)
	{
	  TICK(); waxpby(nrow, 1.0, r, 0.0, r, p); TOCK(t2);
	}
      else
	{
	  oldrtrans = rtrans;
	  TICK(); ddot (nrow, r, r, &rtrans, t4); TOCK(t1);// 2*nrow ops
    AD_INTERMEDIATE(rtrans, "r_dot"); 
	  AD_real beta = rtrans/oldrtrans;
	  TICK(); waxpby (nrow, 1.0, r, beta, p, p);  TOCK(t2);// 2*nrow ops
	  for(int ii=0; ii < nrow; ii++) {AD_INTERMEDIATE(p[ii], "p_wxpy");}; 
	}
      normr = sqrt(rtrans);
      if (rank==0 && (k%print_freq == 0 || k+1 == max_iter))
      cout << "Iteration = "<< k << "   Residual = "<< AD_value(normr) << endl;
     

#ifdef USING_MPI
      TICK(); exchange_externals(A,p); TOCK(t5); 
#endif
      TICK(); HPC_sparsemv(A, p, Ap); TOCK(t3); // 2*nnz ops

	for(int ii=0; ii < nrow; ii++) {AD_INTERMEDIATE(Ap[ii], "matvec");}; 
      AD_real alpha = 0.0;
      TICK(); ddot(nrow, p, Ap, &alpha, t4); TOCK(t1); AD_INTERMEDIATE(alpha, "a_dot");
      alpha = rtrans/alpha;
      TICK(); waxpby(nrow, 1.0, x, alpha, p, x);// 2*nrow ops
	for(int ii=0; ii < nrow; ii++) {AD_INTERMEDIATE(x[ii], "x_wxpy");}; 
      waxpby(nrow, 1.0, r, -alpha, Ap, r);  TOCK(t2);// 2*nrow ops
	for(int ii=0; ii < nrow; ii++) {AD_INTERMEDIATE(r[ii], "r_wxpy");}; 
      niters = k;
    }


  // Store times
  times[1] = t1; // ddot time
  times[2] = t2; // waxpby time
  times[3] = t3; // sparsemv time
  times[4] = t4; // AllReduce time
#ifdef USING_MPI
  times[5] = t5; // exchange boundary time
#endif
  delete [] p;
  delete [] Ap;
  delete [] r;
  times[0] = mytimer() - t_begin;  // Total time. All done...
  return(0);
}
