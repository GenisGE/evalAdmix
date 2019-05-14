#include <math.h>
#include <fstream>

void adaptF(double **f, int **G, int K, double **q, int N, int M, int withoutInd, int nIts, double **aNorm, double **bNorm){


  //  double tol=0.00001;

for(int it=0; it<nIts;it++){



for(int i=0;i<N;i++){
  for(int j=0;j<M;j++){
    double fpart=0;
    for(int k=0;k<K;k++){
      fpart=fpart+f[j][k]*q[i][k];
      //      fprintf(stderr, "added %f times %f to fpart", f[j][k], q[i][k]);
    }
    if(fpart==0){
	fprintf(stderr, "0 in abNorm ind %d site %d",i,j);
	exit(0);
      }
    aNorm[i][j]=1.0/fpart;
    fpart=1-fpart;
    bNorm[i][j]=1.0/fpart;
    //    if(isnan(aNorm[i][j])||isnan(bNorm[i][j])){
    if(aNorm[i][j]==0){
	fprintf(stderr, "0 in abNorm ind %d site %d",i,j);
	exit(0);
      }
  }

}
 
 for(int j=0;j<M;j++){
    for(int k=0;k<K;k++){
      double sumAG=0;
      double sumBG=0;
      for(int i=0;i<N;i++){
        if(i==withoutInd)
           continue;
        sumAG=sumAG+G[i][j]*f[j][k]*q[i][k]*aNorm[i][j];
        sumBG=sumBG+(2-G[i][j])*q[i][k]*(1-f[j][k])*bNorm[i][j];
	if(isnan(sumAG)){
	    fprintf(stderr, "NaN produced when multiplying %d, %f,%f and %f",G[i][j],f[j][k],q[i][k],aNorm[i][j]);
	    exit(0);
}
      }
      f[j][k]=sumAG/(sumAG+sumBG);
      if(isnan(f[j][k])){
	fprintf(stderr, "Nan in f site %d k %d, have sumAG = %f, sumBG = %f", j,k, sumAG, sumBG);
	  exit(0);
	}
    }
  }



 /* 
for(int k=0;k<K;k++){
  for(int j=0;j<M;j++){
    //   f[j*K[0]+k] = fnew[j*K[0]+k];
    if(f[j][k] < tol)
       f[j][k] = tol;
    if(f[j][k] > 1- tol)
       f[j][k] = 1- tol;

  }
}
 */

}

}
