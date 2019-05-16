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


eachPars *allPars = NULL;


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
  fprintf(stderr,"\t-autosomeMax 22\t autosome ends with this chromsome\n");
  fprintf(stderr,"\t-nIts 5\t number of iterations to do for frequency correction\n");



  exit(0);

}



void *functionC(void *a) //the a means nothing
{
  int running_job;

  pthread_mutex_lock(&mutex1); // Protect so only one thread can edit jobs done
  
  while (jobs > 0) {
    running_job = jobs--;
    pthread_mutex_unlock(&mutex1); // Allow threads to go crazy

    ////////////////////////////////////////////// not protected
    int c = NumJobs-running_job;
    eachPars p=allPars[c];
    myPars *pars=p.pars;
    int i=p.ind1;
evalAdmix(p.cor,p.ind1,pars->r,pars->mean_r,pars->data->matrix,pars->K,pars->Q,pars->F,pars->nInd,pars->nSites,pars->nIts, pars->isMissing);
    
    //////////////////////////////////////////////


  }

  return NULL;
}



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
  
  const char *outname = "output.k";
  int autosomeMax = 23;
  int nIts = 5;
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
    else if (strcmp(argv[argPos],"-a")==0){
      autosomeMax = atoi(argv[argPos+1])+1; 
    }
    else if(strcmp(argv[argPos],"-fname")==0 || strcmp(argv[argPos],"-f")==0) 
      fname=argv[argPos+1]; 
    else if(strcmp(argv[argPos],"-qname")==0 || strcmp(argv[argPos],"-q")==0) 
      qname=argv[argPos+1];
    else if(strcmp(argv[argPos],"-nThreads")==0 || strcmp(argv[argPos],"-P")==0) 
      nThreads=atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-nIts")==0)
      nIts=atoi(argv[argPos+1]);
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
  
  //check if files exists
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

  pars->K=K;
  pars->nSites=nSites;
  pars->nInd=nInd;
  
  fp=fopen(outname,"w");
  double **F =allocDouble(nSites,K);
  double **Q =allocDouble(nInd,K);
  pars->F=F;
  pars->Q=Q;
  readDouble(Q,nInd,K,qname,0);
  readDoubleGZ(F,nSites,K,fname,1);

 
    double *cor=new double[(nInd-1)*nInd/2];
    double **pi=allocDouble(nInd, nSites);
    double **r=allocDouble(nInd,nSites);
    double *mean_r=new double[nInd];


    // Create matrix to avoid missing data
    int **isMissing= new int*[nInd];
    for(size_t i=0; i<nInd; i++)
      isMissing[i] = new int[nSites];

    for(int i=0; i < nInd; i++){
      for(int j=0; j< nSites; j++){
	if(pars->data->matrix[i][j]==3){
	  isMissing[i][j] = 1;
	}
	else{
	  isMissing[i][j] = 0;
	}

      }
    }
    
    // Calculate unadapted residuals
    calcRes(pi, r, mean_r, pars->data->matrix, pars-> Q, pars->F, K, nSites, nInd, isMissing);
    fprintf(stderr, "Finished calculating normal residuals\n");


    pars -> r = r;
    pars -> mean_r = mean_r;
    pars -> nIts = nIts;
    pars -> isMissing = isMissing;
    // without threading
    
    if(nThreads==1){
      // Loop over inds, adapt its freq and calculate cor with the remaining inds
      for(int i=0;i<(nInd-1);i++){

	fprintf(stderr, "Estimating freqs without ind %d\r",i);
	
	evalAdmix(cor, i, pars -> r, pars -> mean_r, pars->data->matrix,pars->K,pars-> Q, pars->F,pars->nInd, pars->nSites,pars->nIts, pars->isMissing);

      }

 
    }
    else{ // with threads (the cool way)

      fprintf(stderr, "Correcting frequencies with %d threads...",nThreads);
      
      NumJobs = nInd-1;
      jobs =  nInd-1;
      cunt = 0;
      allPars = new eachPars[NumJobs];

    
      for(int c=0;c<NumJobs;c++){ // fill allPars with each job
    
      
	allPars[c].ind1=c;
	allPars[c].pars=pars;
	allPars[c].cor=cor;
      }
      pthread_t thread1[nThreads];
    
      for (int i = 0; i < nThreads; i++)
	pthread_create(&thread1[i], NULL, &functionC, NULL);
    
      // Wait all threads to finish
      for (int i = 0; i < nThreads; i++)
	pthread_join(thread1[i], NULL);
    
     

    }


    fprintf(stderr, "\nFinished, going to write all correlations\n");
    
    ///////// print header
    fprintf(fp,"ind1\tind2\tcor_res\n");

    int idxWrite = 0;
    for(int i=0;i<(nInd-1);i++){
      for(int j=i+1;j<nInd;j++){
	fprintf(fp,"%d\t%d\t%f\n",i,j,cor[idxWrite]);
	idxWrite ++;
      }

    }


  


// clean


 for(int j = 0; j < nSites; j++){ 
    delete[] F[j];
    }


 for(int i = 0; i < nInd; i++){
   delete[] Q[i];
     //delete[] pi[i];
   delete[] r[i];
 }
 delete[] F;
 delete[] Q;
   //delete[] pi;
 delete[] r;
 delete[] mean_r;
 fclose(fp);
 delete[] allPars;
 

  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fprintf(stderr, "\t[ALL done] results have been outputted to %s\n",outname);



  return(0);
}

