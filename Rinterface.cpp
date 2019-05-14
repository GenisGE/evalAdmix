#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "evalAdmix.h"
#include <stdio.h>
#include <string>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <vector>

#ifndef _types_h
#define _types_h
#include "types.h"
#endif

#ifndef _alloc_h
#define _alloc_h
#include "alloc.h"
#endif

#include "filereader_and_conversions.h"
#include "r_c_conversions.h"
using namespace std;

extern "C" {SEXP relateC (SEXP Rfreq,SEXP Rstart,SEXP Rgeno1,SEXP Rgeno2,SEXP Ra1,SEXP Ra2,SEXP RuseSq,SEXP RtolStop,SEXP RnIter,SEXP RK,SEXP Rtol);
 }


SEXP relateC (SEXP Rfreq,SEXP Rstart,SEXP Rgeno1,SEXP Rgeno2,SEXP Ra1,SEXP Ra2,SEXP RuseSq,SEXP RtolStop,SEXP RnIter,SEXP RK,SEXP Rtol) { 

  //stuff to return
  SEXP ans = R_NilValue;
  SEXP ans_name = R_NilValue;
  PROTECT(ans    = allocVector(VECSXP, 2)); //numer of stuff to return
  PROTECT(ans_name = allocVector(STRSXP, 2));//added 5

  //ugly R to C stuff
  double tolStop;
  int nSites,K,nIter,useSq,numIter; //number of sites
  int *geno1,*geno2;
  double *a1,*a2,*start;
  double **f;
  PROTECT(Rgeno1);
  nSites=length(Rgeno1);
  UNPROTECT(1);
  geno1 = r2cIntArray(Rgeno1);
  geno2 = r2cIntArray(Rgeno2);
  a1 = r2cDoubleArray(Ra1);
  a2 = r2cDoubleArray(Ra2);
  f = r2cDoubleMatrix(Rfreq);
  K=INTEGER(RK)[0];
  nIter=INTEGER(RnIter)[0];
  tolStop=REAL(RtolStop)[0];
  start = r2cDoubleArray(Rstart);
  useSq=INTEGER(RuseSq)[0];
  double tol = REAL(Rtol)[0];

  //run da shit
 evalAdmix(tolStop,nSites,K,nIter,useSq,numIter,geno1,geno2,a1,a2,start,f,tol);


 // clean
 delete[] geno1;
 delete[] geno2;
 for(int i=0;i<nSites;i++)
   delete[] f[i];
 delete[] f;
 delete[] a1;
 delete[] a2;

 //ugly C to R fuck
 SEXP theK,theIter;
 double k0=start[0];
 double k1=start[1];
 double k2=start[2];
 delete[] start;
 PROTECT(theK=allocVector(REALSXP,3));
 REAL(theK)[0]=k0;
 REAL(theK)[1]=k1;
 REAL(theK)[2]=k2;
 UNPROTECT(1);
     
 PROTECT(theIter=allocVector(INTSXP,1));
 INTEGER(theIter)[0]=numIter;
 UNPROTECT(1);
 SET_VECTOR_ELT(ans, 0, theK);
 SET_STRING_ELT(ans_name,0 , mkChar("k"));
 SET_VECTOR_ELT(ans, 1, theIter);
 SET_STRING_ELT(ans_name,1, mkChar("nIter"));
 setAttrib(ans, R_NamesSymbol, ans_name);
 UNPROTECT(2);
 return(ans) ;
} 

