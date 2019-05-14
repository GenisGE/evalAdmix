#include <iostream>
#include <fstream>
#include "Cinterface.h"
#include <string>
#include <iomanip>//used for setting number of decimals in posterior
#include <cstdlib>
#include <zlib.h>
#include <sstream>
#include <cstring>
#include <vector>
#include <sys/stat.h>
#ifndef _types_h
#define _types_h
#include "types.h"
#endif
#include <pthread.h>
#include <assert.h>

#ifndef _alloc_h
#define _alloc_h
#include "alloc.h"
#endif





#include "filereader_and_conversions.h"
#include "extractors.h"
#include "asort.h"
#include "evalAdmix.h"

using namespace std;

pthread_t *threads = NULL;
pthread_t *threads1 = NULL;
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
int NumJobs;
int jobs;
int *printArray;
FILE *fp;
int cunt;
int doInbreeding =0;

eachPars *allPars = NULL;
/*
void *relateWrap(void *a){
  eachPars *p = (eachPars *)a;
  myPars *pars=p->pars;
  int i=p->ind1;
  int j=p->ind2;
  int numIt=0;
  
  fprintf(stderr,"%f\n",pars->tolStop);
  fprintf(stderr,"%f\n",pars->tol);
  fprintf(stderr,"%d\n",pars->nSites);
  fprintf(stderr,"%d\n",i);
  fprintf(stderr,"%d\n",j);
  fprintf(stderr,"%d\n",pars->K);
  fprintf(stderr,"%d\n",pars->maxIter);
  fprintf(stderr,"%d\n",pars->useSq);
  fprintf(stderr,"%d\n",pars->data->matrix[i][0]);
 fprintf(stderr,"%d\n",pars->data->matrix[j][0]);
 fprintf(stderr,"%f\n",pars->Q[j][0]);
fprintf(stderr,"%f\n",p->start[0]);
  
  //if(doInbreeding==0)
  evalAdmix(pars->tolStop,pars->nSites,pars->K,pars->maxIter,pars->useSq,numIt,pars->data->matrix[i],pars->data->matrix[j],pars->Q[i],pars->Q[j],p->start,pars->F,pars->tol);
  //  else
  //    ibAdmix(pars->tolStop,pars->nSites,pars->K,pars->maxIter,pars->useSq,numIt,pars->data,pars->Q,p->start,pars->F,pars->tol,i,pars->likes);
  p->numIter=numIt;
  return NULL;
}
*/
void readDoubleGZ(double **d,int x,int y,const char*fname,int neg){
  fprintf(stderr,"opening : %s with x=%d y=%d\n",fname,x,y);
  const char*delims=" \n";
  gzFile fp = NULL;
  if((fp=gzopen(fname,"r"))==NULL){
    fprintf(stderr,"cant open:%s\n",fname);
    exit(0);
  }
  int lens=1000000 ;
  char buf[lens];
  for(int i=0;i<x;i++){
    if(NULL==gzgets(fp,buf,lens)){
 	fprintf(stderr,"Error: Only %d sites in frequency file (maybe increase buffer)\n",i);
	exit(0);
    }
    if(neg)
      d[i][0] = 1-atof(strtok(buf,delims));
    else
      d[i][0] = atof(strtok(buf,delims));
    for(int j=1;j<y;j++){
      //      fprintf(stderr,"i=%d j=%d\n",i,j);
      if(neg)
	d[i][j] = 1-atof(strtok(NULL,delims));
      else
	d[i][j] = atof(strtok(NULL,delims));
    }
  }
  if(NULL!=gzgets(fp,buf,lens)){
    fprintf(stderr,"Error: Too many sites in frequency file. Only %d sites in plink file\n",x);
    exit(0);

  }
  gzclose(fp);
}

void readDouble(double **d,int x,int y,const char*fname,int neg){
  fprintf(stderr,"opening : %s with x=%d y=%d\n",fname,x,y);
  const char*delims=" \n";
  FILE *fp = NULL;
  if((fp=fopen(fname,"r"))==NULL){
    fprintf(stderr,"can't open:%s\n",fname);
    exit(0);
  }
  int lens=1000000 ;
  char buf[lens];
  for(int i=0;i<x;i++){
    if(NULL==fgets(buf,lens,fp)){
    	fprintf(stderr,"Error: Only %d individuals in admxiture (-qname) file\n",i);
      exit(0);
    }
    if(neg)
      d[i][0] = 1-atof(strtok(buf,delims));
    else
      d[i][0] = atof(strtok(buf,delims));
    for(int j=1;j<y;j++){
      //      fprintf(stderr,"i=%d j=%d\n",i,j);
      if(neg)
	d[i][j] = 1-atof(strtok(NULL,delims));
      else
	d[i][j] = atof(strtok(NULL,delims));
    }
  }
   if(NULL!=fgets(buf,lens,fp)){
    fprintf(stderr,"Error: Too many individuals in admixture (-qname) file. Only %d individuals in plink file\n",x);
    exit(0);
  }

  fclose(fp); 
}

int getK(const char*fname){
  const char*delims=" \n";
  FILE *fp = NULL;
  if((fp=fopen(fname,"r"))==NULL){
    fprintf(stderr,"can't open:%s\n",fname);
    exit(0);
  }
  int lens=100000 ;
  char buf[lens];
  if(NULL==fgets(buf,lens,fp)){
    fprintf(stderr,"Increase buffer\n");
    exit(0);
  }
  strtok(buf,delims);
  int K=1;
  while(strtok(NULL,delims))
    K++;
  fclose(fp);

  return(K);
}


double **allocDouble(size_t x,size_t y){
  double **ret= new double*[x];
  for(size_t i=0;i<x;i++)
    ret[i] = new double[y];
  return ret;
}


double ***allocDouble3D(size_t x,size_t y,size_t z){
  double ***ret= new double**[x];
  for(size_t i=0;i<x;i++){
    ret[i] = new double*[y];
    for(size_t j=0; j<y;j++)
      ret[i][j] = new double[z];
  }
    return ret;
}



void info(){
  fprintf(stderr,"Arguments:\n");
  fprintf(stderr,"\t-plink name of the binary plink file (excluding the .bed)\n");
  fprintf(stderr,"\t-fname Ancestral population frequencies\n"); 
  fprintf(stderr,"\t-qname Admixture proportions\n"); 
  fprintf(stderr,"\t-o name of the output file\n"); 

  fprintf(stderr,"Setup:\n"); 
  fprintf(stderr,"\t-P Number of threads\n");
  fprintf(stderr,"\t-F 1\t if you want to estimate inbreeding\n"); 
  fprintf(stderr,"\t-autosomeMax 22\t autosome ends with this chromsome\n"); 



  exit(0);

}


/*
void *functionC(void *a) //the a means nothing
{
  int running_job;

  pthread_mutex_lock(&mutex1);
  
  while (jobs > 0) {
    running_job = jobs--;
    pthread_mutex_unlock(&mutex1);

    ////////////////////////////////////////////// not protected
    int c = NumJobs-running_job;
    eachPars p=allPars[c];
    myPars *pars=p.pars;
    int i=p.ind1;
    int j=p.ind2;
    int numIt=0;
    evalAdmix(pars->tolStop,pars->nSites,pars->K,pars->maxIter,pars->useSq,numIt,pars->data->matrix[i],pars->data->matrix[j],pars->Q[i],pars->Q[j],p.start,pars->F,pars->tol);
    p.numIter=numIt;
    p.numI[0]=numIt;
    
    //////////////////////////////////////////////

    pthread_mutex_lock(&mutex1);

    int d = NumJobs-running_job;
    printArray[d]=1;
    if(d%50==0)
      fprintf(stderr,"\rrunning i1:%d i2:%d",allPars[d].ind1,allPars[d].ind2);  

    while(cunt<NumJobs){
    
      if(printArray[cunt]==0)
	break;

      fprintf(fp,"%d\t%d\t%f\t%f\t%f\t%d\n",allPars[cunt].ind1,allPars[cunt].ind2,allPars[cunt].start[0],allPars[cunt].start[1],allPars[cunt].start[2],allPars[cunt].numI[0]);
      //fprintf(fp,"%d\t%d\t%f\t%f\t%f\t%d\t%d\n",allPars[cunt].ind1,allPars[cunt].ind2,allPars[cunt].start[0],allPars[cunt].start[1],allPars[cunt].start[2],allPars[cunt].numI[0],cunt);
        cunt++;
    }

  }
  pthread_mutex_unlock(&mutex1);

  return NULL;
}
*/


void fex(const char* fileName){
  const char*delims=" \n";
  FILE *fp = NULL;
  if((fp=fopen(fileName,"r"))==NULL){
    fprintf(stderr,"can't open:%s\n",fileName);
    exit(0);
  }
  fclose(fp);
}


int main(int argc, char *argv[]){
  clock_t t=clock();//how long time does the run take
  time_t t2=time(NULL);
  
  // print commandline
  for(int i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  
  // if not arguments are given print information
  if(argc==1){
    info();
    return 0;
 }
  
  int useSq=1;
  double tolStop=0.000001;
  int maxIter=5000;
  
  int numIter=4000;
  double tol=0.0001;
  const char *outname = "output.k";
  int autosomeMax = 23;
  string geno= "";
  string pos = "";
  string chr = "";
  string plink_fam;
  string plink_bim;
  string plink_bed;
  int nThreads =1;
  bArray *plinkKeep = NULL; //added in 0.97;
  const char* fname = NULL;
  const char* qname = NULL;
  

  //parse arguments
  int argPos=1;
  while(argPos <argc){
    if (strcmp(argv[argPos],"-o")==0){
      outname  = argv[argPos+1]; 
    }
    else if (strcmp(argv[argPos],"-F")==0)
      doInbreeding = atoi(argv[argPos+1]); 
    else if (strcmp(argv[argPos],"-accel")==0)
      useSq = atoi(argv[argPos+1]); 
    else if (strcmp(argv[argPos],"-a")==0){
      autosomeMax = atoi(argv[argPos+1])+1; 
    }
    else if(strcmp(argv[argPos],"-fname")==0 || strcmp(argv[argPos],"-f")==0) 
      fname=argv[argPos+1]; 
    else if(strcmp(argv[argPos],"-qname")==0 || strcmp(argv[argPos],"-q")==0) 
      qname=argv[argPos+1];
    else if(strcmp(argv[argPos],"-nThreads")==0 || strcmp(argv[argPos],"-P")==0) 
      nThreads=atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-plink")==0){
      std::string p_str =string( argv[argPos+1]);
      if(p_str.length()>4){
	std::string ext = p_str.substr(p_str.length()-4,p_str.length());
	if (!ext.compare(".bed")||!ext.compare(".bim")||!ext.compare(".fam")){
	  std::string front = p_str.substr(0,p_str.length()-4);
	  plink_bim = (front+".bim");
	  plink_fam = (front+".fam");
	  plink_bed = (front+".bed");
	}else{
	  plink_bim = (p_str+".bim");
	  plink_fam = (p_str+".fam");
	  plink_bed = (p_str+".bed");	
	}}else{
	plink_bim = (p_str+".bim");
	plink_fam = (p_str+".fam");
	plink_bed = (p_str+".bed");	
     	}
    }
    else{
      printf("\nArgument unknown will exit: %s \n",argv[argPos]);
      info();
      return 0;
    }
    
    argPos+=2;
  }
  
  //check if files exits
  // fexists
  fex(qname);
  fex(fname);
 
  //////////////////////////////////////////////////
  //read plink data
  printf("\t-> Will assume these are the plink files:\n\t\tbed: %s\n\t\tbim: %s\n\t\tfam: %s\n",plink_bed.c_str(),plink_bim.c_str(),plink_fam.c_str());
  int numInds = numberOfLines(plink_fam.c_str())-1;//the number of individuals is just the number of lines

  
  myPars *pars =  new myPars();
  plinkKeep = doBimFile(pars,plink_bim.c_str()," \t",autosomeMax);  
  fprintf(stdout,"\t-> Plink file contains %d autosomale SNPs\n",plinkKeep->numTrue);
  fprintf(stdout,"\t-> reading genotypes ");
  fflush(stdout);
  iMatrix *tmp = bed_to_iMatrix(plink_bed.c_str(),numInds,plinkKeep->x);
  fprintf(stdout," - done \n");
  fflush(stdout);
  if(tmp->y==plinkKeep->numTrue){
    
    pars->data = tmp;
  }
  else{
    fprintf(stdout,"\t-> extractOK (%d %d) ",tmp->x,plinkKeep->numTrue);
    fflush(stdout);
    
    pars->data = extractOK(plinkKeep,tmp);
    killMatrix(tmp);
    fprintf(stdout," - done \n");
    fflush(stdout);

  }
  killArray(plinkKeep);
  fprintf(stdout,"\t-> sorting ");
  fflush(stdout);
  mysort(pars,0);
  fprintf(stdout," - done \n");
  fflush(stdout);
    // printf("Dimension of genodata:=(%d,%d), positions:=%d, chromosomes:=%d\n",pars->data->x,pars->data->y,pars->position->x,pars->chr->x);
  if(pars->data->y != pars->chr->x || pars->position->x != pars->data->y){
    printf("Dimension of data input doesn't have compatible dimensions, program will exit\n");
    printf("Dimension of genodata:=(%d,%d), positions:=%d, chromosomes:=%d\n",pars->data->x,pars->data->y,pars->position->x,pars->chr->x);
    return 0;
  }
  int K=getK(qname);
  int nSites=pars->data->y;
  int nInd=pars->data->x;
  fprintf(stderr,"\t\t->K=%d\tnSites=%d\tnInd=%d\n",K,nSites,nInd);
  pars->maxIter=maxIter;
  pars->tol=tol;
  pars->tolStop=tolStop;
  pars->K=K;
  pars->nSites=nSites;
  pars->useSq=useSq;
  
  fp=fopen(outname,"w");
  double **F =allocDouble(nSites,K);
  double **Q =allocDouble(nInd,K);
  pars->F=F;
  pars->Q=Q;
  readDouble(Q,nInd,K,qname,0);
  readDoubleGZ(F,nSites,K,fname,1);

  ///////// print header
  fprintf(fp,"ind1\tind2\tcor_res\n");

  // without threading

  //  if(nThreads==1){
  double *cor=new double[(nInd-1)*nInd/2];
    double **pi=allocDouble(nInd, nSites);
    double **r=allocDouble(nInd,nSites);
    double *mean_r=new double[nInd];
    int nIts=5;
    /*
    double *sumQ=new double[K]; // TRY PRIMING NEWf ESTIMATION

    for(int k=0; k<K;k++){
      for(int i=0;i<nInd;i++){
	sumQ[k] += Q[i][k];
      }
    }
    */
    fprintf(stderr,"running i1:0 i2:0");  

    
  // Calculate unadapted residuals
  calcRes(pi, r, mean_r, pars->data->matrix, pars-> Q, pars->F, K, nSites, nInd);
  fprintf(stderr, "finished calculating normal residuals\n");
  // Loop over inds, adapt its freq and calculate cor with the remaining inds
  for(int i=0;i<(nInd-1);i++){

    evalAdmix(cor, i, r, mean_r, pars->data->matrix,K,pars-> Q, pars->F,nInd, nSites,nIts,sumQ);

      }
  fprintf(stderr, "going to write all correlations\n");
  //fprintf(sterr, "there are in total %d correlations to write", sizeof(cor))
  int idxWrite = 0;
  for(int i=0;i<(nInd-1);i++){
    for(int j=i+1;j<nInd;j++){
      fprintf(fp,"%d\t%d\t%f\n",i,j,cor[idxWrite]);
      idxWrite ++;
    }

  }
      fprintf(stderr, "correlations have been written\n");


  //}
 
    //  }
  /*  else{ // with threads (the cool way)

    NumJobs = nInd*(nInd-1)/2;
    jobs =  nInd*(nInd-1)/2;
    cunt = 0;
    allPars = new eachPars[NumJobs];

    int *indMatrix = new int[nInd*(nInd-1)];
    int cunter=0;
    for(int i=0;i<nInd-1;i++){
      for(int j=i+1;j<nInd;j++){
	indMatrix[cunter*2]=i;
	indMatrix[cunter*2+1]=j;
	cunter++;
      }
    }

    printArray=new int[NumJobs];
    for(int c=0;c<NumJobs;c++){
      printArray[c]=0;
    
      double *start=new double[3];
      int *numI=new int[1];
      start[0]=0.7;
      start[1]=0.2;
      start[2]=0.1;
      int i=indMatrix[c*2];
      int j=indMatrix[c*2+1];
      
      allPars[c].start=start;
      allPars[c].ind1=i;
      allPars[c].ind2=j;
      allPars[c].numIter=0;
      allPars[c].numI=numI;
      allPars[c].pars=pars;
      
    }
    pthread_t thread1[nThreads];
    
    for (int i = 0; i < nThreads; i++)
      pthread_create(&thread1[i], NULL, &functionC, NULL);
    
  // Wait all threads to finish
    for (int i = 0; i < nThreads; i++)
      pthread_join(thread1[i], NULL);
    
    //  for (int c = 0; c < NumJobs; c++)
  // fprintf(fp,"%d\t%d\t%f\t%f\t%f\t%d\n",allPars[c].ind1,allPars[c].ind2,allPars[c].start[0],allPars[c].start[1],allPars[c].start[2],allPars[c].numI[0]);
  

  delete[] indMatrix;
 }
  */


// clean
   fprintf(stderr, "gonna do final clean-up\n");

  fprintf(stderr,"\n");

 for(int j = 0; j < nSites; j++){ 
    delete[] F[j];
    }

  fprintf(stderr, "cleaned inner F\n");

 for(int i = 0; i < nInd; i++){
   delete[] Q[i];
     fprintf(stderr, "cleant one inner Q\n");
     //delete[] pi[i];
     //fprintf(stderr, "cleant one inner pi\n");
   delete[] r[i];
     fprintf(stderr, "cleant one inner r\n");
 }
 fprintf(stderr, "cleant inner Q pi and r\n");
 delete[] F;
  fprintf(stderr, "cleant F\n");
 delete[] Q;
   fprintf(stderr, "cleant Q\n");
   //delete[] pi;
   //fprintf(stderr, "cleant pi\n");
 delete[] r;
   fprintf(stderr, "cleant r\n");
 delete[] mean_r;
   fprintf(stderr, "cleant mean_r\n");
 fclose(fp);
   fprintf(stderr, "closed fp\n");
 delete[] allPars;
   fprintf(stderr, "cleant allPars\n");
 

  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fprintf(stderr, "\t[ALL done] results have been outputted to %s\n",outname);



  return(0);
}

