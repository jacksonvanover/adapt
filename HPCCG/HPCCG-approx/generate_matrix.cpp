
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

// Routine to read a sparse matrix, right hand side, initial guess, 
// and exact solution (as computed by a direct solver).

/////////////////////////////////////////////////////////////////////////

// nrow - number of rows of matrix (on this processor)

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <cstdlib>
#include <cstdio>
#include <cassert>

#include <ctime> 
#include <fstream>
#include "generate_matrix.hpp"
#include "HBFormat.hpp" 

using namespace std;


void readHBSMF(std::string input_file, HPC_Sparse_Matrix **A, double **x, double **b, double **xexact){
  int *colptr = NULL;
  int indcrd;
  char *indfmt = NULL;
  ifstream input;
  char *key = NULL;
  int khi;
  int klo;
  char *mxtype = NULL;
  int ncol;
  int neltvl;
  int nnzero;
  int nrhs;
  int nrhsix;
  int nrow;
  int ptrcrd;
  char *ptrfmt = NULL;
  int rhscrd;
  char *rhsfmt = NULL;
  char *rhstyp = NULL;
  int *rowind = NULL;
  char *title = NULL;
  int totcrd;
  int valcrd;
  char *valfmt = NULL;
  double *values = NULL;

  *A = new HPC_Sparse_Matrix; // Allocate matrix struct and fill it

  cout << "  Reading the file '" << input_file << "'.\n";

  input.open ( input_file.c_str ( ) );

  if ( !input )
  {
    cout << "  Error opening the file.\n";
    return;
  }

  hb_header_read ( input, &((*A)->title), &key, &totcrd, &ptrcrd, &indcrd, 
      &valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl, &ptrfmt, 
      &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix );

  colptr = new int[ncol+1];

  if ( mxtype[2] == 'A' )
  {
    rowind = new int[nnzero];
  }
  else if ( mxtype[2] == 'E' )
  {
    rowind = new int[neltvl];
  }
  else
  {
    cout << "  Illegal value of MXTYPE character 3.\n";
    return;
  }

  cout << "  Reading the structure.\n";

  hb_structure_read ( input, ncol, mxtype, nnzero, neltvl, 
      ptrcrd, ptrfmt, indcrd, indfmt, colptr, rowind );

  if ( mxtype[2] == 'A' )
  {
    values = new double[nnzero];
  }
  else if ( mxtype[2] == 'E' )
  {
    values =  new double[neltvl];
  }
  else
  {
    cout << "\n";
    cout << "TEST05 - Warning!\n";
    cout << "  Illegal value of MXTYPE character 3 = " << mxtype[2] << "\n";
    return;
  }

  cout << "  Reading the values.\n";

  hb_values_read ( input, valcrd, mxtype, nnzero, neltvl, valfmt, values );

  input.close ( );

  // This is a symmetric matrix, so I can 
  // reverse things and consider
  // the CCS matrix as a CRS one.

  int local_nrow = ncol;

  (*A)->list_of_vals = values; 
  (*A)->list_of_inds = rowind;  

  (*A)->nnz_in_row  = colptr; 
  (*A)->ptr_to_vals_in_row  = new double*[local_nrow];
  (*A)->ptr_to_inds_in_row  = new int   *[local_nrow];

  for ( int i = 0; i < nnzero; i++){
    rowind[i] -= 1;
  }

  double * curvalptr = (*A)->list_of_vals;
  int * curindptr = (*A)->list_of_inds;
  int total_elements = 0;
  for ( int i = 0; i < ncol; i++){
    colptr[i] = (colptr[i+1] - colptr[i]); 
    total_elements = colptr[i];
  }


  *x = new double[local_nrow];
  *b = new double[local_nrow];
  *xexact = new double[local_nrow];

  for ( int i = 0 ; i < nrow; i++){
    (*A)->ptr_to_vals_in_row[i] = curvalptr;
    (*A)->ptr_to_inds_in_row[i] = curindptr;
    curvalptr += colptr[i];
    curindptr += colptr[i];
    (*x)[i] = 1.0;
    (*b)[i] = 0;
    (*xexact)[i] = 0;
  }

  (*A)->start_row = 0; 
  (*A)->stop_row = nrow;
  (*A)->total_nrow = nrow;
  (*A)->total_nnz = nnzero;
  (*A)->local_nrow = nrow;
  (*A)->local_ncol = ncol;
  (*A)->local_nnz = nnzero;
}

void generate_matrix(int nx, int ny, int nz, HPC_Sparse_Matrix **A, double **x, double **b, double **xexact)

{
#ifdef DEBUG
  int debug = 1;
#else
  int debug = 0;
#endif

#ifdef USING_MPI
  int size, rank; // Number of MPI processes, My process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  int size = 1; // Serial case (not using MPI)
  int rank = 0;
#endif

  srand(time(0));

  *A = new HPC_Sparse_Matrix; // Allocate matrix struct and fill it
  (*A)->title = 0;


  // Set this bool to true if you want a 7-pt stencil instead of a 27 pt stencil
  bool use_7pt_stencil = false;

  int local_nrow = nx*ny*nz; // This is the size of our subblock
  assert(local_nrow>0); // Must have something to work with
  int local_nnz = 27*local_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)

  int total_nrow = local_nrow*size; // Total number of grid points in mesh
  long long total_nnz = 27* (long long) total_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)

  int start_row = local_nrow*rank; // Each processor gets a section of a chimney stack domain
  int stop_row = start_row+local_nrow-1;


  // Allocate arrays that are of length local_nrow
  (*A)->nnz_in_row          = new int[local_nrow];
  (*A)->ptr_to_vals_in_row  = new double*[local_nrow];
  (*A)->ptr_to_inds_in_row  = new int   *[local_nrow];
  (*A)->ptr_to_diags        = new double*[local_nrow];

  *x = new double[local_nrow];
  *b = new double[local_nrow];
  *xexact = new double[local_nrow];


  // Allocate arrays that are of length local_nnz
  (*A)->list_of_vals = new double[local_nnz];
  (*A)->list_of_inds = new int   [local_nnz];

  double * curvalptr = (*A)->list_of_vals;
  int * curindptr = (*A)->list_of_inds;

  long long nnzglobal = 0;
  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
        int curlocalrow = iz*nx*ny+iy*nx+ix;
        int currow = start_row+iz*nx*ny+iy*nx+ix;
        int nnzrow = 0;
        (*A)->ptr_to_vals_in_row[curlocalrow] = curvalptr;
        (*A)->ptr_to_inds_in_row[curlocalrow] = curindptr;
        for (int sz=-1; sz<=1; sz++) {
          for (int sy=-1; sy<=1; sy++) {
            for (int sx=-1; sx<=1; sx++) {
              int curcol = currow+sz*nx*ny+sy*nx+sx;
              //            Since we have a stack of nx by ny by nz domains , stacking in the z direction, we check to see
              //            if sx and sy are reaching outside of the domain, while the check for the curcol being valid
              //            is sufficient to check the z values
              if ((ix+sx>=0) && (ix+sx<nx) && (iy+sy>=0) && (iy+sy<ny) && (curcol>=0 && curcol<total_nrow)) {
                if (!use_7pt_stencil || (sz*sz+sy*sy+sx*sx<=1)) { // This logic will skip over point that are not part of a 7-pt stencil
                  if (curcol==currow) {
                    (*A)->ptr_to_diags[curlocalrow] = curvalptr;
                    *curvalptr++ = 27.0;
                  }
                  else {
                    *curvalptr++ = -1.0;
                  }
                  *curindptr++ = curcol;
                  nnzrow++;
                } 
              }
            } // end sx loop
          } // end sy loop
        } // end sz loop
        (*A)->nnz_in_row[curlocalrow] = nnzrow;
        nnzglobal += nnzrow;
        (*x)[curlocalrow] = 0.0;
        (*b)[curlocalrow] = 27.0 - ((double) (nnzrow-1)) + (double)(rand()) / (double)(RAND_MAX);
        (*xexact)[curlocalrow] = 1.0;
      } // end ix loop
    } // end iy loop
  } // end iz loop  
  if (debug) cout << "Process "<<rank<<" of "<<size<<" has "<<local_nrow;

  if (debug) cout << " rows. Global rows "<< start_row
    <<" through "<< stop_row <<endl;

  if (debug) cout << "Process "<<rank<<" of "<<size
    <<" has "<<local_nnz<<" nonzeros."<<endl;

  (*A)->start_row = start_row ; 
  (*A)->stop_row = stop_row;
  (*A)->total_nrow = total_nrow;
  (*A)->total_nnz = total_nnz;
  (*A)->local_nrow = local_nrow;
  (*A)->local_ncol = local_nrow;
  (*A)->local_nnz = local_nnz;

  return;
}
