
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
#include <vector>
#include "mytimer.hpp"
#include "HPCCG.hpp"

#define X_WXPY 0
#define P_WXPY 1
#define R_DOT 2
#define R_WXPY 3
#define MATVEC 4
#define A_DOT 5

#define TICK()  t0 = mytimer() // Use TICK and TOCK to time a code section
#define TOCK(t) t += mytimer() - t0

int HPCCG(HPC_Sparse_Matrix * A,
	  const double * const b, double * const x,
	  const int max_iter, const double tolerance, int &niters, double & normr,
	  double * times, double const * relative_err, std::vector<std::pair<double,double>> residual_bounds)

{
  bool approx_flag;
  double t_begin = mytimer();  // Start timing right away

  double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;
#ifdef USING_MPI
  double t5 = 0.0;
#endif
  int nrow = A->local_nrow;
  int ncol = A->local_ncol;

  double * r = new double [nrow];
  double * p = new double [ncol]; // In parallel case, A is rectangular
  double * Ap = new double [nrow];

  normr = 0.0;
  double rtrans = 0.0;
  double oldrtrans = 0.0;

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
  normr = sqrt(rtrans);
  approx_flag = false;
  for (const auto& b : residual_bounds){
    if (normr > b.first && normr < b.second){
      approx_flag = true;
      break;
    }
  }

  if (rank==0) cout << "Initial Residual = "<< normr << endl;

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
    if (approx_flag){
      if (relative_err[R_DOT] > 0) {
        rtrans = rtrans + relative_err[R_DOT] * rtrans;
      }
      else if (relative_err[R_DOT] < 0){
        rtrans = (float) rtrans;
      }
    }
	  
    double beta = rtrans/oldrtrans;
	  TICK(); waxpby (nrow, 1.0, r, beta, p, p);  TOCK(t2);// 2*nrow ops
    if (approx_flag) {
      if (relative_err[P_WXPY] > 0){
        for(int ii=0; ii < nrow; ii++){
          p[ii] = p[ii] + relative_err[P_WXPY] * p[ii];
        }
      }
      else if (relative_err[P_WXPY] < 0){
        for(int ii=0; ii < nrow; ii++){
          p[ii] = (float) p[ii];
        }
      }
    }
	}
      normr = sqrt(rtrans);
      if (rank==0 && (k%print_freq == 0 || k+1 == max_iter))
      cout << "Iteration = "<< k << "   Residual = "<< normr << endl;
      for (const auto& b : residual_bounds){
        if (normr > b.first && normr < b.second){
          approx_flag = true;
          break;
        }
      }

#ifdef USING_MPI
      TICK(); exchange_externals(A,p); TOCK(t5); 
#endif
      TICK(); HPC_sparsemv(A, p, Ap); TOCK(t3); // 2*nnz ops
      if (approx_flag){
        if (relative_err[MATVEC] > 0) {
          for(int ii=0; ii < nrow; ii++){
            Ap[ii] = Ap[ii] + relative_err[MATVEC] * Ap[ii];
          }
        }
        else if (relative_err[MATVEC] < 0){
          for(int ii=0; ii < nrow; ii++){
            Ap[ii] = (float) Ap[ii];
          }
        }
      }
      double alpha = 0.0;
      TICK(); ddot(nrow, p, Ap, &alpha, t4); TOCK(t1);
      if (approx_flag) {
        if (relative_err[A_DOT] > 0) {
          alpha = alpha + relative_err[A_DOT] * alpha;
        }
        else if (relative_err[A_DOT] < 0){
          alpha = (float) alpha;
        }
      }
      alpha = rtrans/alpha;
      TICK(); waxpby(nrow, 1.0, x, alpha, p, x);// 2*nrow ops
      if (approx_flag){
        if (relative_err[X_WXPY] > 0) {
          for(int ii=0; ii < nrow; ii++){
            x[ii] = x[ii] + relative_err[X_WXPY] * x[ii];
          }
        }
        else if (relative_err[X_WXPY] < 0){
          for(int ii=0; ii < nrow; ii++){
            x[ii] = (float) x[ii];
          }          
        }
      }
      waxpby(nrow, 1.0, r, -alpha, Ap, r);  TOCK(t2);// 2*nrow ops
      if (approx_flag){
        if (relative_err[R_WXPY] > 0) {
          for(int ii=0; ii < nrow; ii++){
            r[ii] = r[ii] + relative_err[R_WXPY] * r[ii];
          }
        }
        else if (relative_err[R_WXPY] < 0){
          for(int ii=0; ii < nrow; ii++){
            r[ii] = (float) r[ii];
          }
        }
      }
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
