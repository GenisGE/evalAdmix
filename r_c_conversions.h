int r_int_to_c_int(SEXP i);
double r_real_to_c_double(SEXP i);
int r_bool_to_c_int(SEXP i);
SEXP c_int_to_r_bool(int i);
SEXP c_int_to_r_int(int i);
SEXP c_double_to_r_real(double i);
dArray *r_array_to_c_dArray(SEXP vec);
iArray *r_array_to_c_iArray(SEXP vec);
SEXP c_iArray_to_r_array(iArray *vec);
SEXP c_dArray_to_r_array(dArray *vec);
iMatrix *r_matrix_to_c_iMatrix(SEXP v);
SEXP c_iMatrix_to_r_matrix(iMatrix *v);
SEXP c_dMatrix_to_r_matrix(dMatrix *v);
double *r2cDoubleArray(SEXP vec);
int *r2cIntArray(SEXP vec);
double **r2cDoubleMatrix(SEXP v);
int **r2cIntMatrix(SEXP v);

