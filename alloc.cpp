/*
  Copyright (C) 2008 Thorfinn Sand Korneliussen  thorfinn@binf.ku.dk
  
  This file is part of HMMld 0.881
  
  HMMld is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Foobar is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with HMMld.  If not, see <http://www.gnu.org/licenses/>.
*/  

#include <iostream>
#include <exception>  
#include <cstdlib>
#ifndef _types_h
#define _types_h
#include "types.h"
#endif



#include "alloc.h"

//used for filling up the internal datastructures with variables before using them
#define _fillup_ 1

iMatrix *allocIntMatrix(int x, int y){
  try{
    iMatrix *tmp = new iMatrix();
    int **ppi = new int*[x];
    

    // int *curPtr = new int [x * y];
    for( int i = 0; i < x; ++i) {
      *(ppi + i) = new int[y];
      // curPtr += y;
    }
#if _fillup_
    for (int i=0;i<x;i++)
      for(int n=0;n<y;n++)
	ppi[i][n]=0;
#endif
    printf("\t-> Allocated: %.3f gig memory\n",(float)sizeof(int)*x*y/1000000000);
    tmp->matrix=ppi;
    tmp->x=x;
    tmp->y=y;
    return tmp;
  }catch(std::exception & e){
    printf("\t-> Tried to allocate: %.3f gig memory\n",(float)sizeof(int)*x*y/1000000000);
    printf("\t-> Allocation failed, likely cause: need more memory, will exit.\n");
    exit(0);
  }
}

usiMatrix *allocUSIntMatrix(int x, int y){
  try{
    usiMatrix *tmp = new usiMatrix();
    unsigned short int **ppi = new unsigned short int*[x];
    // int *curPtr = new int [x * y];
    for( int i = 0; i < x; ++i) {
      *(ppi + i) = new unsigned short int[y];
      // curPtr += y;
    }
#if _fillup_
    for (int i=0;i<x;i++)
      for(int n=0;n<y;n++)
	ppi[i][n]=0;
#endif
    printf("\t-> Allocated: %.3f gig memory\n",(float)sizeof(unsigned short int)*x*y/1000000000);
    tmp->matrix=ppi;
    tmp->x=x;
    tmp->y=y;
    return tmp;
  }catch(std::exception & e){
    printf("\t-> Tried to allocate: %.3f gig memory\n",(float)sizeof(int)*x*y/1000000000);
    printf("\t-> Allocation failed, likely cause: need more memory, will exit.\n");
    exit(0);
  }
}



void killArray(dArray *var){
  delete [] var->array;
  delete var;
}


void killArray(bArray *var){
  delete [] var->array;
  delete var;
}

void killArray(iArray *var){
  delete [] var->array;
  delete var;
}


void killDoubleMatrix(double **var){
  delete [] *var;
  delete [] var;
}


void killIntMatrix(int **var){
  delete [] *var;
  delete [] var;
}

void killUSIntMatrix(unsigned short int **var){
  delete [] *var;
  delete [] var;
}

void killMatrix(dMatrix *var){
  killDoubleMatrix(var->matrix);
  delete var;
}


void killMatrix(iMatrix *var){
  killIntMatrix(var->matrix);
  delete var;
}

void killMatrix(usiMatrix *var){
  killUSIntMatrix(var->matrix);
  delete var;
}

void killSnpMatrix(snpMatrix *mat){
  killMatrix(mat->pba);
  killMatrix(mat->dprime);
  killMatrix(mat->pbA);
  killMatrix(mat->pBa);
  killMatrix(mat->pBA);
  killMatrix(mat->rmisc);
  killMatrix(mat->D);
  killMatrix(mat->lod);
  delete mat;
}



void killPars(pars *mat){
  killArray(mat->Sk1);
  killArray(mat->Sk2);
  killArray(mat->Sk3);
  killArray(mat->t);
  killMatrix(mat->S);
  killArray(mat->maf);
  killArray(mat->maf2);
  killArray(mat->ind1);
  killArray(mat->ind2);
  //killArray(mat->pos);
  //killArray(mat->chr);
  delete mat;
}






dMatrix *allocDoubleMatrix(int x, int y){
  try{
    dMatrix *tmp = new dMatrix();
    double **ppi = new double*[x];
    double *curPtr = new double [x * y];
    
    for( int i = 0; i < x; ++i) {
      *(ppi + i) = curPtr;
      curPtr += y;
    }
#if _fillup_
    for (int i=0;i<x;i++)
      for(int n=0;n<y;n++)
	ppi[i][n]=0;
#endif
    tmp->matrix=ppi;
    tmp->x=x;
    tmp->y=y;
    return tmp;
  }catch(std::exception & e){
        printf("\t-> Tried to allocate: %.3f gig memory\n",(float)sizeof(double)*x*y/1000000000);
    printf("\t-> Allocation failed, likely cause: need more memory, will exit.\n");
    exit(0);
  }
}

/*
  Same as above but with an extra argument
  this will be on every entry of the matrix.

*/
dMatrix *allocDoubleMatrix(int x, int y,float inputVar){
  dMatrix *retVal = allocDoubleMatrix(x,y);
  for (int i=0;i<x;i++)
    for(int n=0;n<y;n++)
      retVal->matrix[i][n]= inputVar;

  return retVal;
}



iArray *allocIntArray(size_t  i){
  try{
    int *tmp= new int[i];
    iArray *intArray = new iArray();
#if _fillup_
    for (unsigned int j=0;j<i;j++)
      tmp[j]=0;
#endif
    intArray->x=i;
    intArray->array=tmp;
    return intArray;
  }catch(std::exception & e){
    printf("\t-> Tried to allocate: %.3f gig memory\n",(float)sizeof(size_t)*i/1000000000);
    printf("\t-> Thats an array of length: %lu \n",i);
    printf("\t-> Allocation failed, likely cause: need more memory, will exit.\n");
    exit(0);
  }
}


bArray *allocBoolArray(int i){
  try{
    int *tmp= new int[i];
  bArray *boolArray = new bArray();
#if _fillup_
  for (int j=0;j<i;j++)
    tmp[j]=0;
#endif
  boolArray->x=i;
  boolArray->array=tmp;
  boolArray->numTrue = 0;
  return boolArray;
 }catch(std::exception & e){
    printf("\t-> Allocation failed, likely cause: need more memory, will exit.\n");
    exit(0);
  }

}


dArray *allocDoubleArray(int i){
  try{
    double *tmp= new double[i];
    dArray *darray = new dArray();
#if _fillup_
    for (int j=0;j<i;j++)
      tmp[j]=0.0;
#endif
    darray->x=i;
    darray->array=tmp;
    return darray;
  }catch(std::exception & e){
    printf("\t-> Tried to allocate: %.3f gig memory\n",(float)sizeof(double)*i/1000000000);
    printf("\t-> Allocation failed, likely cause: need more memory, will exit.\n");
    exit(0);
  }
}
//copy constructor
dArray *allocDoubleArray(dArray *other){
  try{
    dArray *darray = new dArray();
  double *tmp= new double[other->x];
  
  for (int j=0;j<other->x;j++)
    tmp[j]=other->array[j];
  darray->x=other->x;
  darray->array=tmp;
  return darray;
  }catch(std::exception & e){
    printf("\t-> Allocation failed, likely cause: need more memory, will exit.\n");
    exit(0);
  }

}


void killFunctionPars(functionPars *tmp){

  killMatrix(tmp->data);
  killArray(tmp->pair);
  
  killArray(tmp->position);//this is the pos from the outer run
  killArray(tmp->alim);
  killArray(tmp->par); 
  killArray(tmp->chr);

  delete tmp;

}

void killHapStruct(hapStruct *tmp){
  killMatrix(tmp->pBA);
  killMatrix(tmp->pBa);
  killMatrix(tmp->pbA);
  killMatrix(tmp->pba);
  killMatrix(tmp->mea);
  delete tmp;
}

void killFullResults(fullResults *tmp){

  killArray(tmp->kResult);

  if (tmp->S1!=tmp->S)
    killMatrix(tmp->S);
  killMatrix(tmp->S1);

  killArray(tmp->t);

  
  killArray(tmp->position);
  
  //killArray(tmp->alim);
  
  if(tmp->choose!=NULL)
    killArray(tmp->choose);
  if(tmp->back!=0)
    killHapStruct(tmp->hap);
  else
    delete tmp->hap;
  killArray(tmp->maf);
  killArray(tmp->chr); //COKE ROX
  killMatrix(tmp->post);
  if(tmp->convInfo!=NULL)
    killMatrix(tmp->convInfo);
  killArray(tmp->keepList);
  killArray(tmp->path);
  delete tmp;
}


void printFunctionPars(const functionPars *pars){
  printf("----------------------------------------\n");
  if(pars!=NULL){
    printf("\tWill start HMMld with following parameters\n");
    printf("\tdoParameter=%d,\tpair[0]=%d,\tpair[1]=%d,\n",pars->doPar,pars->pair->array[0],pars->pair->array[1]);
    printf("\tdoAllpairs=%d,",pars->doAllPairs);
    if(pars->doPar!=0)
      printf("\tdoParameter calulation with values (%f\t,%f\t,%f)\n",pars->par->array[0],pars->par->array[1],pars->par->array[2]);
    printf("\talim=[%.3f,\t%.3f],\n",pars->alim->array[0],pars->alim->array[1]);
    printf("\tfixK2=%d,\tfixA=%d,\t\tcalcA=%d,\n",pars->fixK,pars->fixA,pars->calcA);
    printf("\tfixK_val=%.3f,\tfixA_value=%.3f,\n",pars->fixK_val,pars->fixA_val);
    printf("\tdoPrune=%d,\tprune_val=%.3f,\n",pars->doPrune,pars->prune_val);
    printf("\tmin=%.3f,\teps=%.3f,\tphi=%.3f,\tconvTol=%.3f\n",pars->min,pars->epsilon,pars->phi,pars->convTol);
    printf("\tdouble_recom=%d\t,LD=%d\t,ld_adj=%d,\tback=%d,\tback2=%d\n",pars->double_recom,pars->LD,pars->ld_adj,pars->back,pars->back2);
    printf("\ttimesToConverge=%d\t,timesToRun=%d\n",pars->timesToConverge,pars->timesToRun);
    //printf("\tdoChromo:%d\tdoMoment:%d\n",pars->doChromo,pars->moment);
    printf("\tpost_file_output:%s \t k_file_output:%s \n",pars->postFilename.c_str(),pars->kFilename.c_str());
    printf("\tWill print results: %d\n",pars->print_results);
  }
  printf("----------------------------------------\n");
}

void flush_print(const char prefix[],const char msg[]){
  printf("%s%s",prefix,msg);
  fflush(stdout);
}

void flush_print(const char prefix[],const char msg[],double val){
  printf("%s%s%f",prefix,msg,val);
  fflush(stdout);
}

void flush_print(const char prefix[],const char msg[],int val){
  printf("%s%s%d",prefix,msg,val);
  fflush(stdout);
}


void flush_print(const char msg[]){
  char pre[] = "\t-> ";
  flush_print(pre,msg);
}

void flush_print(const char msg[],double val){
  char pre[] = "\t-> ";
  flush_print(pre,msg,val);
}

void flush_print(const char msg[],int val){
  char pre[] = "\t-> ";
  flush_print(pre,msg,val);
}

//added in 0.98


fArray *allocFloatArray(int i){
  try{
    float *tmp= new float[i];
    fArray *darray = new fArray();
#if _fillup_
    for (int j=0;j<i;j++)
      tmp[j]=0.0;
#endif
    darray->x=i;
    darray->array=tmp;
    return darray;
  }catch(std::exception & e){
printf("\t-> Tried to allocate: %.3f gig memory\n",(float)sizeof(float)*i/1000000000);
    printf("\t-> Allocation failed, likely cause: need more memory, will exit.\n");
    exit(0);
  }

}


fMatrix *allocFloatMatrix(size_t x, size_t y){
  try{
    fMatrix *tmp = new fMatrix();
    float **ppi = new float*[x];
    float *curPtr = new float [x * y];
    
    for(unsigned long int i = 0; i < x; ++i) {
      *(ppi + i) = curPtr;
      curPtr += y;
    }
#if _fillup_
    for (unsigned long int i=0;i<x;i++)
      for(unsigned long int n=0;n<y;n++)
	ppi[i][n]=0;
#endif
    
    tmp->matrix=ppi;
    tmp->x=x;
    tmp->y=y;
    return tmp;
  }catch(std::exception & e){
    printf("\t-> Tried to allocate: %.3f gig memory\n",(float)sizeof(float)*x*y/1000000000);
    printf("\t-> Allocation failed, likely cause: need more memory, will exit.\n");
    exit(0);
  }

}



void killArray(fArray *var){
  if(var!=NULL){
    delete [] var->array;
    delete var;
  }
}


void killFloatMatrix(float **var){
  if(var!=NULL){
    delete [] *var;
    delete [] var;
  }
}



void killMatrix(fMatrix *var){
  if(var!=NULL){
    killFloatMatrix(var->matrix);
    delete var;
  }
}
