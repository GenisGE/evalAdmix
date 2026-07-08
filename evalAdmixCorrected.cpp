#include "evalAdmixCorrected.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>

// Invert a small K x K matrix in place. Returns true on success.
static bool invertKxK(double *A, int K) {
  double *T = new double[K * 2 * K];
  for (int i = 0; i < K; i++)
    for (int j = 0; j < K; j++) {
      T[i*2*K + j]     = A[i*K + j];
      T[i*2*K + K + j] = (i == j) ? 1.0 : 0.0;
    }
  for (int i = 0; i < K; i++) {
    int piv = i;
    double maxv = std::fabs(T[i*2*K + i]);
    for (int r = i + 1; r < K; r++) {
      double v = std::fabs(T[r*2*K + i]);
      if (v > maxv) { maxv = v; piv = r; }
    }
    if (maxv == 0.0) { delete[] T; return false; }
    if (piv != i)
      for (int c = 0; c < 2*K; c++)
        std::swap(T[i*2*K + c], T[piv*2*K + c]);
    double inv = 1.0 / T[i*2*K + i];
    for (int c = 0; c < 2*K; c++) T[i*2*K + c] *= inv;
    for (int r = 0; r < K; r++) {
      if (r == i) continue;
      double f = T[r*2*K + i];
      if (f == 0.0) continue;
      for (int c = 0; c < 2*K; c++) T[r*2*K + c] -= f * T[i*2*K + c];
    }
  }
  for (int i = 0; i < K; i++)
    for (int j = 0; j < K; j++) A[i*K + j] = T[i*2*K + K + j];
  delete[] T;
  return true;
}

static void mulMatVecKxK(const double *A, const double *v, double *out, int K) {
  for (int i = 0; i < K; i++) {
    double s = 0;
    for (int j = 0; j < K; j++) s += A[i*K + j] * v[j];
    out[i] = s;
  }
}

void evalAdmixCorrected(double **cor,
                        unsigned short int **G,
                        int **isMissing,
                        double **Qmat,
                        int K,
                        int nInd,
                        int nSites,
                        int nIndUse,
                        int *useInds)
{
  // Mapping subset index -> raw individual index
  int *orig = new int[nIndUse];
  for (int i = 0, j = 0; i < nInd; i++)
    if (useInds[i]) orig[j++] = i;

  // Qused = Q restricted to useInds (nIndUse x K)
  double **Qused = new double*[nIndUse];
  for (int i = 0; i < nIndUse; i++) {
    Qused[i] = new double[K];
    for (int k = 0; k < K; k++) Qused[i][k] = Qmat[orig[i]][k];
  }

  // Per-individual mean genotype mu_i and heterozygosity d_i over non-missing sites
  double *muInd = new double[nIndUse];
  double *dInd  = new double[nIndUse];
  int    *nUsed = new int[nIndUse];
  for (int i = 0; i < nIndUse; i++) { muInd[i] = 0; dInd[i] = 0; nUsed[i] = 0; }
  for (int j = 0; j < nSites; j++)
    for (int i = 0; i < nIndUse; i++) {
      int raw = orig[i];
      if (!isMissing[j][raw]) {
        double g = (double)G[j][raw];
        muInd[i] += g;
        dInd[i]  += g * (2.0 - g);
        nUsed[i]++;
      }
    }
  for (int i = 0; i < nIndUse; i++) {
    if (nUsed[i] > 0) { muInd[i] /= nUsed[i]; dInd[i] /= nUsed[i]; }
    else { muInd[i] = 0; dInd[i] = 0; }
  }

  // Q^T Q (K x K)
  double *QtQ = new double[K*K];
  for (int a = 0; a < K; a++)
    for (int b = 0; b < K; b++) {
      double s = 0;
      for (int i = 0; i < nIndUse; i++) s += Qused[i][a] * Qused[i][b];
      QtQ[a*K + b] = s;
    }
  double *QtQinv = new double[K*K];
  std::memcpy(QtQinv, QtQ, sizeof(double)*K*K);
  if (!invertKxK(QtQinv, K)) {
    fprintf(stderr, "Error: Q^T Q is singular; cannot compute corrected correlation.\n");
    exit(1);
  }

  // P (nIndUse x nIndUse) projector. P[i][j] = Qused[i] . (QtQinv Qused[j])
  double **P = new double*[nIndUse];
  for (int i = 0; i < nIndUse; i++) P[i] = new double[nIndUse];
  double *yj = new double[K];
  for (int j = 0; j < nIndUse; j++) {
    mulMatVecKxK(QtQinv, Qused[j], yj, K);
    for (int i = 0; i < nIndUse; i++) {
      double dot = 0;
      for (int k = 0; k < K; k++) dot += Qused[i][k] * yj[k];
      P[i][j] = dot;
    }
  }
  delete[] yj;

  // Residual matrix R [nSites][nIndUse]
  // R[j][i] = g[j][i] - (P g)[j][i]; missing genotypes mean-imputed
  //           only inside the projection, residual set to NaN if missing.
  float  **R = new float*[nSites];
  double *gimp = new double[nIndUse];
  double *yyj  = new double[K];
  double *zcol = new double[K];
  for (int j = 0; j < nSites; j++) {
    R[j] = new float[nIndUse];
    for (int i = 0; i < nIndUse; i++) {
      int raw = orig[i];
      gimp[i] = isMissing[j][raw] ? muInd[i] : (double)G[j][raw];
    }
    for (int k = 0; k < K; k++) {
      double s = 0;
      for (int i = 0; i < nIndUse; i++) s += Qused[i][k] * gimp[i];
      zcol[k] = s;
    }
    mulMatVecKxK(QtQinv, zcol, yyj, K);
    for (int i = 0; i < nIndUse; i++) {
      double pred = 0;
      for (int k = 0; k < K; k++) pred += Qused[i][k] * yyj[k];
      int raw = orig[i];
      R[j][i] = isMissing[j][raw] ? NAN : (float)((double)G[j][raw] - pred);
    }
  }
  delete[] gimp;
  delete[] yyj;
  delete[] zcol;

  // Per-individual mean residual (over non-missing sites)
  double *meanR = new double[nIndUse];
  for (int i = 0; i < nIndUse; i++) {
    double s = 0; int c = 0;
    for (int j = 0; j < nSites; j++)
      if (!std::isnan(R[j][i])) { s += R[j][i]; c++; }
    meanR[i] = (c > 0) ? s / c : 0.0;
  }

  // b_hat: pairwise Pearson correlation with pairwise site skipping,
  // using individual means (matches original correlateRes convention).
  double **bhat = new double*[nIndUse];
  for (int i = 0; i < nIndUse; i++) bhat[i] = new double[nIndUse];
  for (int i = 0; i < nIndUse; i++)
    for (int i2 = i + 1; i2 < nIndUse; i2++) {
      double num = 0, denI = 0, denI2 = 0;
      for (int j = 0; j < nSites; j++) {
        if (!isMissing[j][orig[i]] && !isMissing[j][orig[i2]]) {
          double ri  = (double)R[j][i]  - meanR[i];
          double ri2 = (double)R[j][i2] - meanR[i2];
          num  += ri * ri2;
          denI += ri * ri;
          denI2+= ri2 * ri2;
        }
      }
      double denom = std::sqrt(denI * denI2);
      bhat[i][i2] = (denom > 0) ? num / denom : 0.0;
    }

  // c_hat = cov2cor of C := (I - P) D_hat (I - P)
  // C[i][i2] = d_i [i==i2] - d_i P[i][i2] - P[i][i2] d_i2 + sum_l P[i][l] d_l P[l][i2]
  double **chat = new double*[nIndUse];
  for (int i = 0; i < nIndUse; i++) chat[i] = new double[nIndUse];
  for (int i = 0; i < nIndUse; i++)
    for (int i2 = 0; i2 < nIndUse; i2++) {
      double v = (i == i2) ? dInd[i] : 0.0;
      v -= dInd[i] * P[i][i2];
      v -= P[i][i2] * dInd[i2];
      double sum = 0;
      for (int l = 0; l < nIndUse; l++) sum += P[i][l] * dInd[l] * P[l][i2];
      v += sum;
      chat[i][i2] = v;
    }
  for (int i = 0; i < nIndUse; i++)
    for (int i2 = i + 1; i2 < nIndUse; i2++) {
      double ci  = chat[i][i];
      double cj  = chat[i2][i2];
      double denom = std::sqrt(ci * cj);
      double r = (denom > 0) ? chat[i][i2] / denom : 0.0;
      cor[i][i2] = bhat[i][i2] - r;
    }

  // Cleanup
  for (int i = 0; i < nIndUse; i++) delete[] Qused[i];  delete[] Qused;
  for (int i = 0; i < nIndUse; i++) delete[] P[i];      delete[] P;
  for (int j = 0; j < nSites; j++) delete[] R[j];      delete[] R;
  for (int i = 0; i < nIndUse; i++) delete[] bhat[i];   delete[] bhat;
  for (int i = 0; i < nIndUse; i++) delete[] chat[i];   delete[] chat;
  delete[] orig;
  delete[] muInd;
  delete[] dInd;
  delete[] nUsed;
  delete[] meanR;
  delete[] QtQ;
  delete[] QtQinv;
}