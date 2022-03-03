/*
  Copyright (C) 2008 Thorfinn Sand Korneliussen  thorfinn@binf.ku.dk
  
  This file is part of HMMld 0.83
  
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
  
//alloc by allocIntMatrix
//deall by killMatrix
//inplmementation of member methods in alloc.cpp
#include <string>
#include <stdio.h>

typedef struct  {
  int x;
  int y;
  int** matrix;
  void info(const std::string name){
    printf("iMatrix:%s dims= (%d,%d) \n",name.c_str(),x,y);
  }

  void print(const char *st,const char *file){
    if(file==NULL){
      if(st!=NULL)
	printf("\niMatrix dims: %s = (%d,%d)\n",st,x,y);
      for (int i=0;i<x;i++){
	for (int j=0;j<y;j++)
	  printf("%d ",matrix[i][j] );
	printf("\n");
      }
    } else{
      FILE *pFile;
      pFile = fopen(file,"w");
      if(st!=NULL)

	fprintf(pFile,"\niMatrix dims of: %s = (%d,%d)\n",st,x,y);
      for (int i=0;i<x;i++){
	for (int j=0;j<y;j++)
	  fprintf(pFile,"%d ",matrix[i][j]);
	fprintf(pFile,"\n");
      }
      fclose(pFile);
    }
  }
}iMatrix ;



typedef struct  {
  int x;
  int y;
  unsigned short int** matrix;
  void info(const std::string name){
    printf("iMatrix:%s dims= (%d,%d) \n",name.c_str(),x,y);
  }

  void print(const char *st,const char *file){
    if(file==NULL){
      if(st!=NULL)
	printf("\niMatrix dims: %s = (%d,%d)\n",st,x,y);
      for (int i=0;i<x;i++){
	for (int j=0;j<y;j++)
	  printf("%hu ",matrix[i][j] );
	printf("\n");
      }
    } else{
      FILE *pFile;
      pFile = fopen(file,"w");
      if(st!=NULL)

	fprintf(pFile,"\niMatrix dims of: %s = (%d,%d)\n",st,x,y);
      for (int i=0;i<x;i++){
	for (int j=0;j<y;j++)
	  fprintf(pFile,"%hu ",matrix[i][j]);
	fprintf(pFile,"\n");
      }
      fclose(pFile);
    }
  }
}usiMatrix ;




//alloc by allocDoubleMatrix
//deall by killMatrix
typedef struct{
  int x;
  int y;
  double** matrix;
  void info(const std::string name) const{
    printf("dMatrix ims :%s = (%d,%d) \n",name.c_str(),x,y);
  }

  void print(const char *st,const char *file){
    if(file==NULL){
      if(st!=NULL)
	printf("\ndMatrix dims of: %s = (%d,%d)\n",st,x,y);
      for (int i=0;i<x;i++){
	for (int j=0;j<y;j++)
	  printf("%f ",matrix[i][j] );//coke_rox
	printf("\n");
      }
    } else{
      FILE *pFile;
      pFile = fopen(file,"w");
      if(st!=NULL)
	fprintf(pFile,"\ndMarix dims of: %s = (%d,%d)\n",st,x,y);
      for (int i=0;i<x;i++){
	for (int j=0;j<y;j++)
	  fprintf(pFile,"%f ",matrix[i][j]);
	fprintf(pFile,"\n");
      }
      fclose(pFile);
    }
  }
} dMatrix ;

//alloc by allocDoubleMatrix
//deall by killMatrix
typedef struct{
  int x;
  int * array;
  void info(const char *name){
    printf("\niArray dims :%s = (%d) \n",name,x);
  }
  void print(const char *st,const char *file){
    if(file==NULL){
      if(st!=NULL)
	printf("\niArray dims: %s = (%d)\n",st,x);
      for (int i=0;i<x;i++)
	printf("%d ",array[i] );
      printf("\n");
      
    } else{
      FILE *pFile;
      pFile = fopen(file,"w");
      if(st!=NULL)
	fprintf(pFile,"\niArray dims of: %s = (%d)\n",st,x);
      for (int i=0;i<x;i++)
	fprintf(pFile,"%d ",array[i]);
      fprintf(pFile,"\n");
      
      fclose(pFile);
    }
  }
}iArray;

typedef struct{
  int x; //length
  int numTrue; //number of Trues;
  int * array;
 void info(const char *name)const{
   printf("bArray dims:%s = (%d) \tnumTrue=%d\n",name,x,numTrue);
  }
  void print(const char *st,const char *file){
    if(file==NULL){
      if(st!=NULL)
	printf("\nbArray dims of: %s = (%d)\tnumTrues=%d\n",st,x,numTrue);
      for (int i=0;i<x;i++)
	printf("%d ",array[i] );
      printf("\n");
      
    } else{
      FILE *pFile;
      pFile = fopen(file,"w");
      if(st!=NULL)
	fprintf(pFile,"\nbArray dims of: %s = (%d)\tnumTrue=%d\n",st,x,numTrue);
      for (int i=0;i<x;i++)
	fprintf(pFile,"%d ",array[i]);
      fprintf(pFile,"\n");
      
      fclose(pFile);
    }
  }
}bArray;

//alloc by allocDoubleArray
//deall by killArray
typedef struct {
  int x;
  double * array;
  void info(const char *name){
    printf("\ndArray dims:%s = (%d) \n",name,x);
  }
  void print(const char *st,const char *file){
    if(file==NULL){
      if(st!=NULL)
	printf("\ndArray dims: %s = (%d)\n",st,x);
      for (int i=0;i<x;i++)
	printf("%f ",array[i] );
      printf("\n");
      
    } else{
      FILE *pFile;
      pFile = fopen(file,"w");
      if(st!=NULL)
	fprintf(pFile,"\ndArray dims: %s = (%d)\n",st,x);
      for (int i=0;i<x;i++)
	fprintf(pFile,"%f ",array[i]);
      fprintf(pFile,"\n");
      
      fclose(pFile);
    }
  }
}dArray;


typedef struct {
  /*
    matrix[d][s]
    d := depth
    s := snp
  */
  dMatrix *dprime;
  dMatrix *pba;
  dMatrix *pbA;
  dMatrix *pBa;
  dMatrix *pBA;
  dMatrix *rmisc;
  dMatrix *D;
  dMatrix *lod;
}snpMatrix;


typedef struct {
  dMatrix *pBA;
  dMatrix *pBa;
  dMatrix *pbA;
  dMatrix *pba;
  dMatrix *mea;
}hapStruct;



typedef struct {
  dArray *Sk1;
  dArray *Sk2;
  dArray *Sk3;
  dArray *t;
  dMatrix *S;
  dMatrix *S1;
  dArray *maf;
  dArray *maf2;
  iArray *ind1;
  iArray *ind2;
  //iArray *chr;
  dArray *pos;
  iArray *choose;
  bArray *usedSnps;
  hapStruct *hap;
  iArray *keepList;
}pars;




typedef struct  { //pars used by relateHMM
  //?? min;
  //pointer
  iMatrix *data; //1
  iArray *pair; //2
  iArray *chr; //1
  dArray *position;//1
  dArray *par; //1
  dArray *alim; //2
  std::string postFilename;
  std::string kFilename;
  // branching vars
  int doAllPairs;
  int doPar; //1
 
 
  int doPrune;
  int back2;
  int fixA;
  int fixK;
  int calcA;//
  int double_recom;
  int LD;
  int ld_adj;
  int print_results;
  // values
  double phi;//
  double calcA_val;
  double min;
  double epsilon;
  int back;
  int timesToRun;
  int timesToConverge;
  double prune_val;
  double fixA_val;
  double fixK_val;
  double convTol;
  int moment;

  std::string stripped_geno;
  std::string stripped_pos;
  std::string stripped_chr;
  std::string stripped_keeplist;

}functionPars ;//pars used by relateHMM


typedef struct {
  dArray *kResult;//0
  double kLike;
  dMatrix *S;
  double kr;
  double a;
  double uLike ; //5
  //conv
  int LD;
  int timesRun;
  int timesConverged;
  dArray *t; //7
  int snp;
  dArray *position;
  int double_recom;//10
  dArray *alim;
  double poLike; //12
  iArray *choose;
  int back;
  hapStruct *hap;
  dMatrix *S1;
  dArray *maf;
  iArray *chr;  //19
  dMatrix *post;
  dMatrix *convInfo;
  bArray *usedSnps;
  iArray *keepList;
  iArray *path;
  dMatrix *mea_global;
  dMatrix *pba_global;
  dMatrix *pBa_global;
  dMatrix *pbA_global;
  dMatrix *pBA_global;
}fullResults;

typedef struct{
  iMatrix *data;
  iArray *chr;
  dArray *pos;
}fromSnpMatrix;

//the following was add in 0.98 from HMMtest 0.73

typedef struct {
  int x;
  float * array;
  void info(const std::string name){
    printf("fArray:%s dims= (%d) \n",name.c_str(),x);
  }
  void print(const char *st,const char *file){
    if(file==NULL){
      printf("\nfArray dims: %s = (%d)\n",st,x);
      for (int i=0;i<x;i++)
	printf("%f ",array[i] );
      printf("\n");
      
    } else{
      FILE *pFile;
      pFile = fopen(file,"w");
      fprintf(pFile,"\nfArray dims of: %s = (%d)\n",st,x);
      for (int i=0;i<x;i++)
	fprintf(pFile,"%f ",array[i]);
      fprintf(pFile,"\n");
      
      fclose(pFile);
    }
  }
}fArray;



typedef struct{
  int x;
  int y;
  float** matrix;
} fMatrix ; 




typedef struct{ //used by HMMtest
  char* postFileName;
  std::string kFileName;
  std::string cFileName;
  std::string xFileName;
  char* dumpPostFileName;
  char* dFileName;
  char* posFileName;
  std::string outFileName;
  int design;
  int missingExists;
  int ccAll;
  int nSim;
  int isK0;
  float sig;
  int infer;
  int every; //use every snp default 1;
}functionPars2; //used by HMMtest


/*
typedef struct  { //pars used by relateHMM
  //?? min;
  //pointer
  iMatrix *data; //1
  iArray *chr; //1
  dArray *position;//1

  double **F;
  double **Q;
  double tol;
  double tolStop;
  int maxIter;
  int useSq;
  int K;
  int nSites;
  double likes;

}myPars ;//pars used by relateHMM
*/



typedef struct  { //pars used by evalAdmix
  //?? min;
  //pointer
  usiMatrix *data; //1
  iArray *chr; //1
  dArray *position;//1
  //int **genos;  transposed geno matrix for faster indexing in loops
  unsigned short int **genos;
  
  double **F;
  double **Q;
  float **r;
  float *mean_r;
  int nIts;
  int K;
  int nInd;
  int nIndUse;
  int nSites;
  double likes;
  int **isMissing;
  int *useInds;

}myPars ;//pars used by evalAdmix


typedef struct  { //pars used by evalAdmix
  //?? min;
  //pointer
  double **genos; //1

  double **F;
  double **Q;
  double **r;
  double *mean_r;
  int nIts;
  int K;
  int nInd;
  int nIndUse;
  int nSites;
  double likes;
  char **keeps;
  int *useInds;

}myNGSPars ;//pars used by evalAdmix ngs version


/*
typedef struct  { //pars used by relateHMM
  //?? min;
  //pointer
  myPars *pars; //1
  int ind1;
  int ind2;
  double *start;
  int numIter;
  int *numI;
  double llh;
}eachPars ;//pars used by relateHMM
 
*/


typedef struct  { //pars used by each job in evalAdmix
  //?? min;
  //pointer
  myPars *pars; //1
  int ind1;
  int ind1u;
  double *cor;
}eachPars ;//pars used by each job in evalAdmix


typedef struct  { //pars used by each job in evalAdmix
  //?? min;
  //pointer
  myNGSPars *pars; //1
  int ind1;
  int ind1u;
  double *cor;
}eachNGSPars ;//pars used by each job in evalAdmix



//some struct will all the data from the beagle file (copied from ngsadmix32.cpp)
typedef struct{
  double **genos; // genotype likelihoods [nSites, 3nInd]
  char *major;
  char *minor;
  char **ids;
  int nSites;
  int nInd;
  char **keeps; //matrix containing 0/1 indicating if data or missing
  int *keepInd; //keepInd[nSites] this is the number if informative samples
  float *mafs;
}bgl;
