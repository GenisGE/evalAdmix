#include <math.h>
#include <fstream>


double **allocDouble2(size_t x,size_t y){
  double **ret= new double*[x];
  for(size_t i=0;i<x;i++)
    ret[i] = new double[y];
  return ret;
}


void adaptF(double **F, int  **G, int K, double **Q, int nInd, int nSites, int withoutInd, int nIts, int **isMissing){


  for(int it=0; it<nIts; it++){
    
    for(int j=0; j<nSites; j++){

	double sumAG[K];
	double sumBG[K];
    
	for(int k=0; k<K;k++){
	  sumAG[k] = 1;
	  sumBG[k] = 1;
	}

	for(int i=0; i<nInd;i++){
	
	  if(i == withoutInd)
	    continue;


	  if(!isMissing[j][i]){
	    double fpart=0;
	    double aNorm;
	    double bNorm;

	    for(int k=0;k<K;k++){
	      fpart+=F[j][k]*Q[i][k];
	    }
	   
	    
	    aNorm = 1.0/fpart;
	    bNorm = 1.0/(1-fpart);

	    for(int k=0; k<K; k++){
	    sumAG[k] += G[j][i]*F[j][k]*Q[i][k]*aNorm; // Minor allele
	    sumBG[k] += (2-G[j][i])*Q[i][k]*(1-F[j][k])*bNorm; // Major allele
	    }
	   
	
	  }
	}

	for(int k=0;k<K;k++)
	  F[j][k] = sumAG[k]/(sumAG[k]+sumBG[k]);
	
    }
  }
}


void calcRes(double **r, double *mean_r, int **G, double **Q, double **F, int K, int nSites, int nInd,int nIndUse,  int **isMissing, int *useInd){

  double sum_r[nIndUse];
  int usedSites[nIndUse];


  for(int i=0; i<nIndUse;i++){
    sum_r[i] = 0;
    usedSites[i] = 0;
  }

  for(int j=0; j<nSites;j++){
    int u=0;
    for(int i=0; i<nInd;i++){
      if(useInd[i]){
      if(!isMissing[j][i]){
	
	double fpart=0; // individual freq pi
	
	for(int k=0; k<K; k++){
	  fpart += Q[i][k]*F[j][k];

	  
	}

	 
	r[u][j] = G[j][i] - 2 * fpart;
	sum_r[u] += r[u][j];
	usedSites[u]++;
	u++;
      }
      }
      
    }
  }

  for(int i=0; i<nIndUse;i++)
    mean_r[i] = sum_r[i]/usedSites[i];

}



 
double correlateRes(int nSites, double *r1, double *r2, double mean_r1, double mean_r2, int **isMissing, int ind1, int ind2){

  double cor_num=0;
  double cor_den1=0;
  double cor_den2=0;
  double d_r1=0;
  double d_r2=0;
  

  for(int j=0; j<nSites;j++){
    if(j==368)
      continue;
    if(!isMissing[j][ind1]&!isMissing[j][ind2]){
      d_r1 = (r1[j] - mean_r1);
      d_r2 = (r2[j] - mean_r2);
      cor_num += d_r1*d_r2;
      cor_den1 += d_r1*d_r1;
      cor_den2 += d_r2*d_r2;
    }
    
  }
  double cor;
  cor = cor_num/(sqrt(cor_den1)*sqrt(cor_den2));
  
  return cor;

}

 

void evalAdmix(double *cor, int ind1u, int ind1, double **r, double *mean_r, int **G, int K, double **Q, double **F, int nInd, int nSites,int nIts, int **isMissing, int nIndUse, int *useInd){

  double **rAdapted=allocDouble2(nIndUse, nSites);
  double *mean_rAdapted = new double[nIndUse];
  double **fAdapted=allocDouble2(nSites,K);


   
  for(int j=0; j<nSites;j++){ // initial values for adapted F are normal F
    for(int k=0; k<K;k++){
      fAdapted[j][k] = F[j][k];
    }
  }

    
  // Estimate new f and calculate new residuals
  adaptF(fAdapted, G, K, Q, nInd, nSites, ind1, nIts,isMissing);
  calcRes(rAdapted, mean_rAdapted, G, Q, fAdapted, K, nSites, nInd, nIndUse, isMissing, useInd);

  for(int ind2=(ind1u+1); ind2<nIndUse; ind2++){
    cor[ind2] =  correlateRes(nSites, r[ind1u], rAdapted[ind2], mean_r[ind1u], mean_rAdapted[ind2], isMissing, ind1, ind2);
  }
  
  // Clean

  for(int j = 0; j < nSites; j++){ 
    delete[] fAdapted[j];

  }

  for(int i = 0; i < nIndUse; i++){
    delete[] rAdapted[i];

  }

 
  delete[] fAdapted;
  delete[] rAdapted;
  delete[] mean_rAdapted;


}
