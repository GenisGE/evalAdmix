
/*
  extractOK and extract
  simply takes column or elements out
 */
//these use a which list
iMatrix *extractOK(iArray *okList,iMatrix *matr);
usiMatrix *extractOK(iArray *okList,usiMatrix *matr);
dMatrix *extractOK(iArray *okList,dMatrix *matr);
iArray *extractOK(iArray *okList,iArray *array);
dArray *extractOK(iArray *okList,dArray *array);

//these use a keeplist
dArray *extractOK(bArray *okList,dArray *array);
iArray *extractOK(bArray *okList,iArray *array);
iMatrix *extractOK(bArray *okList,iMatrix *matr);
usiMatrix *extractOK(bArray *okList,usiMatrix *matr);

//these use a choose list and works only on matrix's
//iArray *extract(iMatrix *var,iArray *choose);
dArray *extract(dMatrix *var,iArray *choose);
iArray *extract(iArray *var,iArray *choose);
