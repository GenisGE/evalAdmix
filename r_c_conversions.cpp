#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#ifndef _types_h
#define _types_h
#include "types.h"
#endif

#ifndef _alloc_h
#define _alloc_h
#include "alloc.h"
#endif


#include "r_c_conversions.h"

int r_int_to_c_int(SEXP i){
  int retVal;
  PROTECT(i = coerceVector(i, INTSXP));
  if(TYPEOF(i) != INTSXP)
    Rprintf("%s expected input is int, you gave wrong type\n",__FUNCTION__);
  retVal = INTEGER(i)[0];
  UNPROTECT(1);
  return retVal;
}

double r_real_to_c_double(SEXP i){
  double retVal;
  PROTECT(i = coerceVector(i, REALSXP));
  if(TYPEOF(i) != REALSXP)
      Rprintf("%s expected input is int, you gave wrong type\n",__FUNCTION__);
  retVal = REAL(i)[0];
  UNPROTECT(1);
  return retVal;
}


int r_bool_to_c_int(SEXP i){
  int retVal =0;
  PROTECT(i = coerceVector(i, LGLSXP));
  if(TYPEOF(i) != LGLSXP)
    Rprintf("%s expected input is bool, you gave wrong type\n",__FUNCTION__);
  if(LOGICAL(i)[0])
    retVal = 1;
  UNPROTECT(1);
  return retVal;
}


SEXP c_int_to_r_bool(int i){
  SEXP retVal = R_NilValue;
  PROTECT(retVal =allocVector(LGLSXP,i));
  if (i!=0)
    LOGICAL(retVal)[0]=TRUE;
  else
    LOGICAL(retVal)[0]=FALSE;
  UNPROTECT(1);
  return retVal;
}

SEXP c_int_to_r_int(int i){
  SEXP retVal = R_NilValue;
  PROTECT(retVal =allocVector(INTSXP,1));
  INTEGER(retVal)[0]=i;

  UNPROTECT(1);
  return retVal;
}

SEXP c_double_to_r_real(double i){
  SEXP retVal = R_NilValue;
  PROTECT(retVal =allocVector(REALSXP,1));
  REAL(retVal)[0] =i;
  UNPROTECT(1);
  return retVal;
}



dArray *r_array_to_c_dArray(SEXP vec){
  PROTECT(vec);
  int len = length(vec);
  dArray *retVal = allocDoubleArray(len);

  for(int i=0;i<len;i++)
    retVal->array[i] = REAL(vec)[i];
  UNPROTECT(1);
  return retVal;
}

iArray *r_array_to_c_iArray(SEXP vec){
  PROTECT(vec);
  int len = length(vec);
  //  printf("len of vec is:%d\n",len);
  iArray *retVal = allocIntArray(len);
  
  for(int i=0;i<len;i++)
    retVal->array[i] = INTEGER(vec)[i];
  UNPROTECT(1);
  return retVal;
}

SEXP c_iArray_to_r_array(iArray *vec){
  SEXP retVal= R_NilValue;
  if(vec==NULL){
    return retVal;
  }
  int len = vec->x;

  PROTECT(retVal=allocVector(INTSXP,len));
  for(int i=0;i<len;i++)
    INTEGER(retVal)[i] = vec->array[i];
  UNPROTECT(1);
  return retVal;
}


SEXP c_dArray_to_r_array(dArray *vec){
  SEXP retVal= R_NilValue;
  if(vec==NULL){
    return retVal;
  }


  int len = vec->x;

  PROTECT(retVal=allocVector(REALSXP,len));
  for(int i=0;i<len;i++){
    REAL(retVal)[i] = vec->array[i];
  }
  UNPROTECT(1);
  return retVal;
}



iMatrix *r_matrix_to_c_iMatrix(SEXP v){
  if(TYPEOF(v) != INTSXP)
    Rprintf("%s input v wrong type\n",__FUNCTION__);
  
  int rows,cols;
  SEXP dim = R_NilValue;
  PROTECT(v);
  PROTECT(dim = getAttrib(v, R_DimSymbol));
  if (length(dim) == 2){
    rows = INTEGER(dim)[0];
    cols = INTEGER(dim)[1];
    //Rprintf("Information: The input contains %i samples with %i snps\n", rows, cols);
  } else
    Rprintf("%s wrong size of dimension of matrix \n",__FUNCTION__);
  iMatrix *retVal = allocIntMatrix(rows,cols);
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++)
      retVal->matrix[i][j] =INTEGER(v)[j*rows+i];
  
  UNPROTECT(2);
  return retVal;
}


SEXP c_iMatrix_to_r_matrix(iMatrix *v){
  SEXP retVal = R_NilValue;
  if(v==NULL){
    return retVal;
  }

  int rows = v->x;
  int cols = v->y;
  // v->info("in sexp");
  PROTECT(retVal=allocMatrix(INTSXP,rows,cols));
  
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++)
      INTEGER(retVal)[j*rows+i] = v->matrix[i][j];
  
  UNPROTECT(1);
  return retVal;
}


SEXP c_dMatrix_to_r_matrix(dMatrix *v){
  SEXP retVal = R_NilValue;
  if(v==NULL){
    return retVal;
  }

  int rows = v->x;
  int cols = v->y;
 
  // v->info("int dmatrix hmmld");
  PROTECT(retVal=allocMatrix(REALSXP,rows,cols));
  
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++){
      REAL(retVal)[j*rows+i] = v->matrix[i][j];
      //      printf("conving %f\t",v->matrix[i][j]);
    }
  UNPROTECT(1);
  return retVal;
}



double *r2cDoubleArray(SEXP vec){
  PROTECT(vec);
  int len = length(vec);
  double *retVal = new double[len];
  
  for(int i=0;i<len;i++)
    retVal[i] = REAL(vec)[i];
  UNPROTECT(1);
  return retVal;
}



int *r2cIntArray(SEXP vec){
  PROTECT(vec);
  int len = length(vec);
  //  printf("len of vec is:%d\n",len);
  int *retVal = new int[len];
  
  for(int i=0;i<len;i++)
    retVal[i] = INTEGER(vec)[i];
  UNPROTECT(1);
  return retVal;
}

int **r2cIntMatrix(SEXP v){
  if(TYPEOF(v) != INTSXP)
    Rprintf("%s input v wrong type\n",__FUNCTION__);
  
  int rows,cols;
  SEXP dim = R_NilValue;
  PROTECT(v);
  PROTECT(dim = getAttrib(v, R_DimSymbol));
  if (length(dim) == 2){
    rows = INTEGER(dim)[0];
    cols = INTEGER(dim)[1];
    //Rprintf("Information: The input contains %i samples with %i snps\n", rows, cols);
  } else
    Rprintf("%s wrong size of dimension of matrix \n",__FUNCTION__);
  int **retVal = new int*[rows];
  for(int i=0;i<rows;i++)
    retVal[i] = new int[cols];

  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++)
      retVal[i][j] =INTEGER(v)[j*rows+i];
  
  UNPROTECT(2);
  return retVal;
}
double **r2cDoubleMatrix(SEXP v){
  if(TYPEOF(v) != REALSXP)
    Rprintf("%s input v wrong type\n",__FUNCTION__);
  
  int rows,cols;
  SEXP dim = R_NilValue;
  PROTECT(v);
  PROTECT(dim = getAttrib(v, R_DimSymbol));
  if (length(dim) == 2){
    rows = INTEGER(dim)[0];
    cols = INTEGER(dim)[1];
    //Rprintf("Information: The input contains %i samples with %i snps\n", rows, cols);
  } else
    Rprintf("%s wrong size of dimension of matrix \n",__FUNCTION__);
  double **retVal = new double*[rows];
  for(int i=0;i<rows;i++)
    retVal[i] = new double[cols];

  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++)
      retVal[i][j] =REAL(v)[j*rows+i];
  
  UNPROTECT(2);
  return retVal;
}
