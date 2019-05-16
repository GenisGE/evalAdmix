#include <math.h>
#include <fstream>


double **allocDouble2(size_t x,size_t y){
  double **ret= new double*[x];
  for(size_t i=0;i<x;i++)
    ret[i] = new double[y];
  return ret;
}


void adaptF(double **F, double **newF, int  **G, int K, double **Q, int N, int M, int withoutInd, int nIts, double **aNorm, double **bNorm, int **isMissing){


  
  for(int it=0; it<nIts;it++){ // do nIts steps to estimate adapted F

    for(int i=0;i<N;i++){ // Calculate normalizing factor for minor (aNorm) major (bNorm) ancestral probability per ind i and site j
      for(int j=0;j<M;j++){
	if(!isMissing[i][j]){
	  double fpart=0;
	  for(int k=0;k<K;k++){
	    fpart=fpart+newF[j][k]*Q[i][k];
	  }

	  if(fpart==0){ // catch some issues I had
	    fprintf(stderr, "0 in abNorm ind %d site %d",i,j);
	    exit(0);
	  }
	  aNorm[i][j]=1.0/fpart;
	  fpart=1-fpart;
	  bNorm[i][j]=1.0/fpart;
	  //    if(isnan(aNorm[i][j])||isnan(bNorm[i][j])){
	  if(aNorm[i][j]==0){ // catch some issues I had
	    fprintf(stderr, "0 in abNorm ind %d site %d",i,j);
	    exit(0);
	  }
	}
      }

    }

    for(int j=0;j<M;j++){
      for(int k=0;k<K;k++){
	double sumAG=0;
	double sumBG=0;
	for(int i=0;i<N;i++){ // Calculate total number of minor (sumAG) and major (sumBG) alleles in each ancestral population
	  if(i==withoutInd)
	    continue;
	  else if(!isMissing[i][j]){
	    sumAG=sumAG+G[i][j]*newF[j][k]*Q[i][k]*aNorm[i][j]; // Minor allele
	    sumBG=sumBG+(2-G[i][j])*Q[i][k]*(1-newF[j][k])*bNorm[i][j]; // Major allele
	    if(isnan(sumAG)){ // catch some other issues i had
	      fprintf(stderr, "NaN produced when multiplying %d, %f,%f and %f",G[i][j],newF[j][k],Q[i][k],aNorm[i][j]);
	      exit(0);
	    }
	  }
	}
	newF[j][k]=sumAG/(sumAG+sumBG); // Calculate new F given new estimate of alleles per K
	if(isnan(newF[j][k])){ // catch more issues i had
	  fprintf(stderr, "Nan in f site %d k %d, have sumAG = %f, sumBG = %f", j,k, sumAG, sumBG);
	  exit(0);
	}
      }


    }
  }
}

void calcRes(double **pi,double **r, double *mean_r, int **g, double **q, double **f, int K, int nSites, int nInd, int **isMissing){

  double sum_r=0;
  
  for(int i=0; i<nInd; i++){
    sum_r = 0;

    for(int j=0; j<nSites;j++){

      if(!isMissing[i][j]){
	pi[i][j]=0;
	for(int k=0;k<K;k++){
	  pi[i][j] += q[i][k]*f[j][k];
	}
	r[i][j] = g[i][j] - 2*pi[i][j];
	sum_r += r[i][j];
      }

    }
    mean_r[i] = sum_r/nSites;
  }
  //fprintf(stderr, "Residuals calculated\n");
}

 
double correlateRes(int nSites, double *r1, double *r2, double mean_r1, double mean_r2, int *isMissing1, int *isMissing2){

  double cor_num=0;
  double cor_den1=0;
  double cor_den2=0;
  double d_r1=0;
  double d_r2=0;
  

  for(int j=0; j<nSites;j++){

    if(!isMissing1[j]&!isMissing2[j]){
      d_r1 = (r1[j] - mean_r1);
      d_r2 = (r2[j] - mean_r2);
      cor_num += d_r1*d_r2;
      cor_den1 += d_r1*d_r1;
      cor_den2 += d_r2*d_r2;
    }
    
  }
  double cor;
  cor = cor_num/(sqrt(cor_den1)*sqrt(cor_den2));
  
  //  fprintf(stderr, "Finished calculating pis without error");
  return cor;

}

 

void evalAdmix(double *cor, int ind1, double **r, double *mean_r, int **G, int K, double **Q, double **F, int nInd, int nSites,int nIts, int **isMissing){

  double **piAdapted=allocDouble2(nInd, nSites);
  double **rAdapted=allocDouble2(nInd, nSites);
  double *mean_rAdapted = new double[nInd];
  double **fAdapted=allocDouble2(nSites,K);
  double **bNorm=allocDouble2(nInd,nSites);
  double **aNorm=allocDouble2(nInd,nSites);


   
  for(int j=0; j<nSites;j++){ // initial values for adapted F are normal F
    for(int k=0; k<K;k++){
      fAdapted[j][k] = F[j][k];
      //newF[j][k] = (F[j][k]*2*init_Nk[k]-G[withoutInd][j]*Q[withoutInd][k])/(2*init_Nk[k]-Q[withoutInd][k]);
    }
  }

    
  // Estimate new f and calculate new residuals
  adaptF(F, fAdapted, G, K, Q, nInd, nSites, ind1, nIts, aNorm, bNorm, isMissing);
  calcRes(piAdapted, rAdapted, mean_rAdapted, G, Q, fAdapted, K, nSites, nInd, isMissing);

  // Index from which to start writing correlation
  int start=ind1*(nInd-1)-(ind1-1)*ind1/2;

  for(int ind2=(ind1+1); ind2<nInd; ind2++){ // Calculate pairwise correlations
    cor[start] = correlateRes(nSites, r[ind1], rAdapted[ind2], mean_r[ind1], mean_rAdapted[ind2], isMissing[ind1], isMissing[ind2]);
    start++;      

  }

  // Clean

  for(int j = 0; j < nSites; j++){ 
    delete[] fAdapted[j];

  }

  for(int i = 0; i < nInd; i++){
    delete[] piAdapted[i];
    delete[] rAdapted[i];
    delete[] aNorm[i];
    delete[] bNorm[i];

  }

 
  delete[] piAdapted;
  delete[] fAdapted;
  delete[] rAdapted;
  delete[] mean_rAdapted;
  delete[] aNorm;
  delete[] bNorm;


}
