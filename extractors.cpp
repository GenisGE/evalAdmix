//contains all the type definition
#ifndef _types_h
#define _types_h
#include "types.h"
#endif

//contains all the allocations
#ifndef _alloc_h
#define _alloc_h
#include "alloc.h"
#endif






iMatrix *extractOK(iArray *okList,iMatrix *matr){
  iMatrix *returnMatrix = allocIntMatrix(matr->x,okList->x);
  for(int j=0;j<matr->x;j++)
    for(int i=0;i<okList->x;i++)
      returnMatrix->matrix[j][i] = matr->matrix [j] [okList->array[i]];
  return returnMatrix;

}

usiMatrix *extractOK(iArray *okList,usiMatrix *matr){
  usiMatrix *returnMatrix = allocUSIntMatrix(matr->x,okList->x);
  for(int j=0;j<matr->x;j++)
    for(int i=0;i<okList->x;i++)
      returnMatrix->matrix[j][i] = matr->matrix [j] [okList->array[i]];
  return returnMatrix;

}


dMatrix *extractOK(iArray *okList,dMatrix *matr){
  dMatrix *returnMatrix = allocDoubleMatrix(matr->x,okList->x);
  for(int j=0;j<matr->x;j++)
    for(int i=0;i<okList->x;i++)
      returnMatrix->matrix[j][i] = matr->matrix [j] [okList->array[i]];
  return returnMatrix;
}

iArray *extractOK(iArray *okList,iArray *array){
  iArray *returnArray = allocIntArray(okList->x);
  for(int j=0;j<okList->x;j++)
    returnArray->array[j] = array->array[okList->array[j]];
  return returnArray;
}

dArray *extractOK(iArray *okList,dArray *array){
  dArray *returnArray = allocDoubleArray(okList->x);
  int extractIndex;
  double var;
  for(int j=0;j<okList->x;j++){
    extractIndex = okList->array[j];
    var = array->array[extractIndex];
    returnArray->array[j] = var;
  }
  return returnArray;
}

dArray *extractOK(bArray *okList,dArray *array){
  dArray *returnArray = allocDoubleArray(okList->numTrue);
  int inPlace = 0;
  for(int j=0 ; j<okList->x ; j++){
    if(okList->array[j] == 1){
      returnArray->array[inPlace] = array->array[j];
      inPlace++;
    }
  }
  return returnArray;
}

iArray *extractOK(bArray *okList,iArray *array){
  iArray *returnArray = allocIntArray(okList->numTrue);
  int inPlace = 0;
  for(int j=0 ; j<okList->x ; j++){
    if(okList->array[j] == 1){
      returnArray->array[inPlace] = array->array[j];
      inPlace++;
    }
  }
  return returnArray;
}

iMatrix *extractOK(bArray *okList,iMatrix *matr){
  iMatrix *returnMatrix = allocIntMatrix(matr->x,okList->numTrue);
  int atPos =0;
  for(int i=0;i<matr->y;i++)
    if(okList->array[i]){
      for(int j=0;j<matr->x;j++)
	returnMatrix->matrix[j][atPos] = matr->matrix [j] [i];
      atPos++;
    }
  return returnMatrix;
}

usiMatrix *extractOK(bArray *okList,usiMatrix *matr){
  usiMatrix *returnMatrix = allocUSIntMatrix(matr->x,okList->numTrue);
  int atPos =0;
  for(int i=0;i<matr->y;i++)
    if(okList->array[i]){
      for(int j=0;j<matr->x;j++)
	returnMatrix->matrix[j][atPos] = matr->matrix [j] [i];
      atPos++;
    }
  return returnMatrix;
}



dArray *extract( dMatrix *var,iArray *choose){
  dArray *retVal = allocDoubleArray(choose->x-1);
   for (int i=1;i<var->y;i++){
    retVal->array[i-1]=var->matrix[choose->array[i]][i];
  }
  return retVal;
}

iArray *extract(iArray *var,iArray *choose){
  iArray *retVal = allocIntArray(var->x-1);
  for (int i=1;i<var->x;i++)
    retVal->array[i-1]=var->array[i-1-choose->array[i]];
  return retVal;
}
