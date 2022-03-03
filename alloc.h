/*
  Copyright (C) 2008 Thorfinn Sand Korneliussen  thorfinn@binf.ku.dk
  
  This file is part of HMMld 0.81
  
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

#ifndef _types_h
#define _types_h
#include "types.h"
#endif


iMatrix *allocIntMatrix(int x, int y);
usiMatrix *allocUSIntMatrix(int x, int y);
void killArray(dArray *var);
void killArray(iArray *var);
void killArray(bArray *var);
void killDoubleMatrix(double **var);
void killIntMatrix(int **var);
void killUSIntMatrix(int **var);
void killMatrix(dMatrix *var);
void killMatrix(iMatrix *var);
void killMatrix(usiMatrix *var);
void killSnpMatrix(snpMatrix *mat);
void killPars(pars *mat);
dMatrix *allocDoubleMatrix(int x, int y);
dMatrix *allocDoubleMatrix(int x, int y,float inputVar);
iArray *allocIntArray(size_t i);
bArray *allocBoolArray(int i);
dArray *allocDoubleArray(int i);
void killFunctionPars(functionPars *tmp);
void killFullResults(fullResults *tmp);
dArray *allocDoubleArray(dArray *other);
void printFunctionPars(const functionPars *pars); //screen output function
void killHapStruct(hapStruct *tmp);


void flush_print(const char prefix[],const char msg[]);
void flush_print(const char msg[]);
void flush_print(const char msg[],double doubleval);
void flush_print(const char msg[],int val);


//added in 0.98
void killMatrix(fMatrix *var);
void killArray(fArray *var);
fMatrix *allocFloatMatrix(size_t x,size_t y);
fArray *allocFloatArray(int i);
