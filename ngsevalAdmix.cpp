#include <math.h>
#include <fstream>
#include <algorithm>

double **allocDouble2b(size_t x,size_t y){
  double **ret= new double*[x];
  for(size_t i=0;i<x;i++)
    ret[i] = new double[y];
  return ret;
}


void adaptFandeG(double **F, double **Q, double **genos, int nSites, int nInd, int K, int withoutInd, char **keep, int nIts){
  
  for(int it=0; it<nIts; it++){
    
    for(int j=0; j<nSites; j++){

      
	
	double sumAG[K];
	double sumBG[K];
    
	for(int k=0; k<K;k++){
	  sumAG[k] = 0;
	  sumBG[k] = 0;
	}

	for(int i=0; i<nInd;i++){
	
	  if(i == withoutInd)
	    continue;

	  if(keep[j][i] == 1){
	    double fpart=0;
	    double aNorm;
	    double bNorm;
	
	    for(int k=0;k<K;k++){
	      fpart+=F[j][k]*Q[i][k];
	    }
	
	    aNorm = 1.0/fpart;
	    bNorm = 1.0/(1-fpart);

	    double pp0=(1-fpart)*(1-fpart)*        genos[j][3*i+2];
	    double pp1=2*(1-fpart)*fpart*  genos[j][3*i+1];
	    double pp2=fpart*fpart*genos[j][3*i+0];
	    double sum=pp0+pp1+pp2;
	    double eG = (pp1+2*pp2)/sum;
	  
	    for(int k=0; k<K;k++){
	      sumAG[k] += eG*F[j][k]*Q[i][k]*aNorm;
	      sumBG[k] += (2-eG)*(1-F[j][k])*Q[i][k]*bNorm;
	    }
	
	  }

	}

	for(int k=0;k<K;k++){
	  F[j][k] = sumAG[k]/(sumAG[k]+sumBG[k]);
	  F[j][k] = std::max(1e-6, std::min(F[j][k],0.999999)); // make sure value is in (1e6, 0.999999) interval 
	    //if(F[j][k]<=0)
	    //F[j][k] = 1e-6;
	}
	
      }    
    }  
}


void ngscalcRes(double **r, double *mean_r, double **F, double **Q, int K, int nSites, int nInd, double **genos, char **keep){

  
  double sum_r[nInd];
  int usedSites[nInd];
  
  for(int i=0; i<nInd;i++){
    sum_r[i] = 0;
    usedSites[i] = 0;
  }
  
  for(int j=0; j<nSites;j++){
    for(int i=0; i<nInd;i++){
       if(keep[j][i] == 1){
	 
	usedSites[i]++;
	
	double fpart=0;
	
	for(int k=0;k<K;k++){
	  fpart+=F[j][k]*Q[i][k];
	  
	} 

	double pp0=(1-fpart)*(1-fpart)*genos[j][3*i+2];
	double pp1=2*(1-fpart)*fpart*genos[j][3*i+1];
	double pp2=fpart*fpart*genos[j][3*i+0];
	double sum=pp0+pp1+pp2;
	double eG=(pp1+2*pp2)/sum;
	
	r[i][j] = eG - 2*fpart;
	/*	if(isnan(r[i][j])){

	    fprintf(stderr, "nan in residual in site %i of individual %i when substracting ind freq of %.3f from expected geno %.3f \n",j,i,eG, fpart);
	    // fprintf(stderr, "fpart is %.3f, sum ");
	    exit(0);
	    
	    }*/
	sum_r[i] += r[i][j];

       }
    }
  }


  for(int i=0; i<nInd;i++)
      mean_r[i] = sum_r[i] / usedSites[i];
}


double correlateNGSRes(int nSites, double *r1, double *r2, double mean_r1, double mean_r2, char **keep, int ind1, int ind2){

  double cor_num=0;
  double cor_den1=0;
  double cor_den2=0;
  double d_r1=0;
  double d_r2=0;
  

  for(int j=0; j<nSites;j++){

    if(keep[j][ind1] == 1 & keep[j][ind2] == 1){
      d_r1 = (r1[j] - mean_r1);
      d_r2 = (r2[j] - mean_r2);

      cor_num += d_r1*d_r2;
      cor_den1 += d_r1*d_r1;
      cor_den2 += d_r2*d_r2;
    }
    
  }
  double cor;
  cor = cor_num/(sqrt(cor_den1)*sqrt(cor_den2));
  /*if(isnan(cor)){
    fprintf(stderr, "NaN generated when calculating residuals of individuals %i and %i \n", ind1, ind2);
    fprintf(stderr, "Numerator was %.3f, denominator was square root of %.3f times %.3f mean residual 1 was %.3f and mean residual 2 %.3f\n",cor_num, cor_den1, cor_den2, mean_r1, mean_r2);
    exit(0);
    }*/
  return cor;

}


void ngsevalAdmix(double *cor, int ind1u, int ind1, double **r, double *mean_r, double **genos, int K, double **Q, double **F, int nInd, int nSites, int nIts, char **keep, int nIndUse){


  double **rAdapted=allocDouble2b(nIndUse, nSites);
  double *mean_rAdapted = new double[nIndUse];
  double **fAdapted=allocDouble2b(nSites,K);


  for(int j=0; j<nSites;j++){ // initial values for adapted F are normal F
    for(int k=0; k<K;k++){
      fAdapted[j][k] = F[j][k];
      //newF[j][k] = (F[j][k]*2*init_Nk[k]-G[withoutInd][j]*Q[withoutInd][k])/(2*init_Nk[k]-Q[withoutInd][k]);
    }
  }

  if(nIts==0){
    
  for(int ind2=(ind1+1); ind2<nInd; ind2++){
    cor[ind2] =  correlateNGSRes(nSites, r[ind1u], r[ind2], mean_r[ind1u], mean_r[ind2], keep, ind1, ind2);
  }

  }
  else{
  // Estimate new f and calculate new residuals
  adaptFandeG(fAdapted, Q, genos, nSites, nInd,K, ind1, keep, nIts);
  ngscalcRes(rAdapted, mean_rAdapted, fAdapted, Q, K, nSites, nIndUse,  genos, keep);
 
  for(int ind2=(ind1+1); ind2<nInd; ind2++){
    cor[ind2] =  correlateNGSRes(nSites, r[ind1], rAdapted[ind2], mean_r[ind1], mean_rAdapted[ind2], keep, ind1, ind2);
  }
  }

  // Clean

  for(int j = 0; j < nSites; j++){ 
    delete[] fAdapted[j];

  }

  for(int i = 0; i < nInd; i++){
    delete[] rAdapted[i];
  }
 
  delete[] fAdapted;
  delete[] rAdapted;
  delete[] mean_rAdapted;
 
}
