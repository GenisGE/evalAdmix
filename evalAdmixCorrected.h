#ifndef _evalAdmixCorrected_h
#define _evalAdmixCorrected_h

// Corrected correlation of residuals (van Waaij et al.) for the admixture model.
// Only supports called genotypes (plink input).
//
// cor       : output [nIndUse][nIndUse]; fills cor[i][j] for j>i (upper triangle).
//             Diagonal is left untouched (writer prints "NA").
// G         : site-major genotypes [nSites][nInd]; values 0/1/2 and 3 = missing.
// isMissing : [nSites][nInd]; 1 if missing (must match G==3).
// Qmat      : admixture proportions [nInd][K].
void evalAdmixCorrected(double **cor,
                        unsigned short int **G,
                        int **isMissing,
                        double **Qmat,
                        int K,
                        int nInd,
                        int nSites,
                        int nIndUse,
                        int *useInds);

#endif