#include <cmath>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "gen_lat_func.h"
extern "C" {
#include "FL.h"
}


/*************************************************************************/
/*************************************************************************/
//GENERAL FUNCTIONS.
/*************************************************************************/
/*************************************************************************/

/*************************************************************************/
//The function CopyAtoB() copies matrix A to matrix B.
//Arguments:
//pA - Pointer to matrix A.
//pB - Pointer to matrix B.
//nel - Total number of elements.
/*************************************************************************/
void CopyAtoB(double* pA, double* pB, int nel) {

  for (int i = 0; i < nel; i++) {
    pB[i] = pA[i];
  }

}

/*************************************************************************/
//The function SqColSums() takes a matrix of reals and returns the column
//sums of the squared matrix. 
//Arguments:
//pmat - Pointer to a matrix.
//pmatSqCS - Pointer to the array of column sums.
//nr - Number of rows of matrix.
//nc - Number of cols of matrix.
/*************************************************************************/
void SqColSums(double* pmat, double* pmatSqCS, int nr, int nc) {

  for (int j = 0; j < nc; j++) {
    pmatSqCS[j] = 0.0;
    for (int i = 0; i < nr; i++) {
      pmatSqCS[j] += pmat[i+j*nr]*pmat[i+j*nr];
    }
  }

}

/*************************************************************************/
//The function SqRowSums() takes a matrix of reals and returns the row
//sums of the squared matrix.
//Arguments:
//pmat - Pointer to a matrix.
//pmatSqRS - Pointer to the array of row sums.
//nr - Number of rows of matrix.
//nc - Number of cols of matrix.
/*************************************************************************/
void SqRowSums(double* pmat, double* pmatSqRS, int nr, int nc) {

  for (int i = 0; i < nr; i++) {
    pmatSqRS[i] = 0.0;
    for (int j = 0; j < nc; j++) {
      pmatSqRS[i] += pmat[i+j*nr]*pmat[i+j*nr];
    }
  }

}

/*************************************************************************/
//The function SqTotSum() takes a matrix of reals and returns the total
//sum of the elements of the squared matrix. 
//Arguments:
//pmat - Pointer to a matrix.
//nel - Total number of elements.
/*************************************************************************/
double SqTotSum(double* pmat, int nel) {

  double ans = 0.0;
  for (int i = 0; i < nel; i++) {
    ans += pmat[i]*pmat[i];
  }
  return(ans);

}

/*************************************************************************/
//The function MatErr() takes pointers to two matrices, poldMat and
//pnewMat, and calculates the square root of the Frobeniuus norm of the
//difference divided by the square root of the Frobenius norm of oldMat.
//Arguments:
//pnewMat - Pointer to the new matrix.
//poldMat - Pointer to the old matrix.
//nel - Total number of elements.
//thresh - The threshold for convergence.  If oldMat=0 and newMat=0,
//function returns 0.  If only oldMat=0, function returns thresh+1.
/*************************************************************************/
double MatErr(double* pnewMat, double* poldMat, int nel, double thresh) {

  double res, oldSqS = SqTotSum(poldMat,nel),
    newSqS = SqTotSum(pnewMat,nel);

  if ((oldSqS==0.0)&&(newSqS==0.0)) {
    res = 0.0;
  } else if (oldSqS==0.0) {
    res = thresh+1.0;
  } else {
    double den = 0.0, num = 0.0;
    for (int i = 0; i < nel; i++) {
      num += (pnewMat[i]-poldMat[i])*(pnewMat[i]-poldMat[i]);
      den += poldMat[i]*poldMat[i];
    }
    res = sqrt(num/den);
  }

  return(res);

}

/*************************************************************************/
//The fuction SoftThresh() applies soft threshholds val by thresh.
//Arguments:
//val - What we want to soft threshold.
//thresh - Thet value we are soft thresholding by.
/*************************************************************************/
double SoftThresh(double val, double thresh) {
  double t1;
  t1 = sign(val);
  double t2;
  t2 = ((fabs(val) - thresh) >= 0) ? fabs(val) - thresh : 0.0;
  return(t1*t2);
}


/*************************************************************************/
/*************************************************************************/
//FUNCTIONS SPECIFIC TO LATENT FEATURE MODEL.
/*************************************************************************/
/*************************************************************************/

/*************************************************************************/
//The function MakeTldY() creates the the vector of \tilde{Y}_j.
//Arguments:
//pnewY - Poimter to array of \tilde{Y}_j.
//pY - Pointer to the original data matrix.
//pB - Pointer to the current Beta matrix.
//pnewT - Pointer to the current Theta matrix, to be updated.
//j - The current row index of Theta.
//S - Number of columns of Y.
//L - Number of rows of Y.
//J - Number of columns of Beta.
/*************************************************************************/
void MakeTldY(double* ptldY, double* pY, double* pB, double* pnewT,
	      int j, int S, int L, int J) {
  for (int s = 0; s < S; s++) {
    ptldY[s] = 0.0;
    for (int l = 0; l < L; l++) {
      double term = BTljsSum(pB,pnewT,s,l,j,S,L,J);
      ptldY[s] += pB[l+j*L]*(pY[l+s*L] - term);
    }
  }
}

/*************************************************************************/
//The function MakeGrvY() creates the the vector of \grave{Y}_j.
//Arguments:
//pgrvY - Pointer to array of \grave{Y}_j.
//pY - Pointer to the original data matrix.
//pT - Pointer to the current Theta matrix.
//pnewB - Pointer to the current Beta matrix, to be updated.
//pTSqRS - Pointer to the array of row sums of Theta^2.
//j - The current row index of Theta.
//S - Number of columns of Y.
//L - Number of rows of Y.
//J - Number of columns of Beta.
/*************************************************************************/
void MakeGrvY(double* pgrvY, double* pY, double* pT, double* pnewB,
	      double* pTSqRS, int j, int S, int L, int J) {
  for (int l = 0; l < L; l++) {
    double term1 = 0.0;
    for (int s = 0; s < S; s++) {
      double term2 = BTljsSum(pnewB,pT,s,l,j,S,L,J);
      term1 += pT[j+s*J]*(pY[l+s*L] - term2);
    }
    pgrvY[l] = term1/pTSqRS[j];
  }
}

/*************************************************************************/
//The function BTljsSum() calculates \sum_{k\ne j}\beta_{lk}\theta_{ks}.
//Arguments:
//pB - Pointer to Beta.
//pT - Pointer to Theta.
//s - Column of Theta.
//l - Row of Beta.
//j - Row of Theta (or column of Beta) not included in sum.
//S - Number of columns of Y.
//L - Number of rows of Y.
//J - Number of columns of Beta.
/*************************************************************************/
double BTljsSum(double* pB, double* pT, int s, int l, int j, int S, int L,
		int J) {
  double res = 0.0;
  for (int k = 0; k < j; k++) {
    res += pB[l+k*L]*pT[k+s*J];
  }
  for (int k = j + 1; k < J; k++) {
    res += pB[l+k*L]*pT[k+s*J];
  }
  return(res);
}

/*************************************************************************/
//The function BC() estimates the value of the Beta matrix for a given
//Theta matrix.  Arguments:
//pnewB - The current Beta which will become updated.
//pY - Pointer to the data matrix.
//pT - Pointer to the current Theta matrix.
//rlam1 - Lambda_1.
//rlam2 - Lambda_2.
//rthresh - The error threshold for determining convergence.
//imaxiter - The maximum number of iterations.
//S - Number of samples.
//L - Number of chromosomal locations.
//J - Number of latent features.
/*************************************************************************/
void BC(double* pnewB, double* pY, double* pT, double rlam1, double rlam2,
	double rthresh, int imaxiter, int S, int L, int J) {

  //The previous Beta to keep track of changes through iterations.
  double *poldB = new double[L*J];
  CopyAtoB(pnewB,poldB,L*J);

  double *pTSqRS = new double[J];
  SqRowSums(pT,pTSqRS,J,S);

  double err = rthresh+1.0, oldBSqS = SqTotSum(poldB,L*J);
  int niter = 0;
  while ((err>rthresh)&&(niter<imaxiter)&&(oldBSqS!=0.0)) {
    UpdateBC(pnewB,pY,pT,pTSqRS,rlam1,rlam2,S,L,J);
    err = MatErr(pnewB,poldB,L*J,rthresh);
    CopyAtoB(pnewB,poldB,L*J);
    oldBSqS = SqTotSum(poldB,L*J);
    niter++;
  }

  delete [] pTSqRS;
  delete [] poldB;

}

/************************************************************************/
//The function UpdateBC() updates Beta.
//Arguments:
//pnewB - Pointer to the current Beta which will be updated.
//pY - Pointer to the original data matrix.
//pT - Pointer to the current Theta matrix.
//pTSqRS - A pointer to the array of row sums of Theta^2.
//rlam1 - Lambda_1.
//rlam2 - Lambda_2.
//S - Number of samples.
//L - Number of chromosomal locations.
//J - Number of latent features.
/************************************************************************/
void UpdateBC(double* pnewB, double* pY, double* pT, double* pTSqRS,
	      double rlam1, double rlam2, int S, int L, int J) {

  int nProt = 0;

  SEXP grvY = R_NilValue, nlam2 = R_NilValue, maxSegs = R_NilValue;
  PROTECT(grvY = allocVector(REALSXP,L));
  nProt++;
  PROTECT(nlam2 = allocVector(REALSXP,1));
  nProt++;
  PROTECT(maxSegs = allocVector(INTSXP,1));
  nProt++;
  INTEGER(maxSegs)[0] = 1000;
  double *pgrvY = REAL(grvY), *pnlam2 = REAL(nlam2);

  for (int j = 0; j < J; j++) {

    if (pTSqRS[j]==0.0) {
      for(int l = 0; l < L; l++) {
	pnewB[l+j*L] = 0.0;
      }
    } else {
      pnlam2[0] = rlam2/pTSqRS[j];
      MakeGrvY(pgrvY,pY,pT,pnewB,pTSqRS,j,S,L,J);
      SEXP retPath = R_NilValue, segmentVec = R_NilValue,
	noNeed = R_NilValue;
      PROTECT(retPath = allocVector(VECSXP,1));
      PROTECT(segmentVec = allocVector(VECSXP,1));
      PROTECT(noNeed =
	      L2L1VitPath(grvY,nlam2,retPath,maxSegs,segmentVec));
      double *pretPath = REAL(VECTOR_ELT(retPath,0));
      int *psegmentVec = INTEGER(VECTOR_ELT(segmentVec,0));
      int nSegs = ncols(VECTOR_ELT(segmentVec,0));
      for (int seg = 0; seg < nSegs; seg++) {
	for (int l = psegmentVec[2*seg]-1; l < psegmentVec[2*seg+1]; l++) {
	  pnewB[l+j*L] = SoftThresh(pretPath[seg],rlam1/(pTSqRS[j]*2.0));
	}
      }
      UNPROTECT(3);
    }
  }

  UNPROTECT(nProt);

}

/***********************************************************************/
//The function TLatL2C estimates the value of the Theta matrix for a
//given Beta matrix.  Arguments:
//pnewT - The current Theta which will be updated. 
//pY - The data matrix.
//pB - The current Beta matrix.
//rthresh - The error threshold for determining convergence.
//imaxiter - The maximum number of iterations.
//rsT -  The constraint on the L2 sum of each row of Theta.
//S - Number of samples.
//L - Number of chromosomal locations.
//J - Number of latent features.
/***********************************************************************/
void TLatL2C(double* pnewT, double* pY, double* pB, double rthresh,
	     int imaxiter, double rsT, int S, int L, int J) {

  //The previous Theta to keep track of changes through iterations.
  double *poldT = new double[J*S];
  CopyAtoB(pnewT,poldT,J*S);

  double *pBSqCS = new double[J];
  SqColSums(pB,pBSqCS,L,J);

  double err = rthresh+1.0;
  int niter = 0;
  while ((err>rthresh)&&(niter<imaxiter)) {
    UpdateTLatL2C(pnewT,pY,pB,pBSqCS,rsT,S,L,J);
    err = MatErr(pnewT,poldT,J*S,rthresh);
    CopyAtoB(pnewT,poldT,J*S);
    niter++;
  }

  delete [] pBSqCS;
  delete [] poldT;

}

/***********************************************************************/
//The function UpdateTLatL2C() updates Theta.
//Arguments:
//pnewT - The matrix which will become the updated Theta.
//pY - Pointer to the original data matrix.
//pB - Pointer to the current Beta matrix.
//pBSqCS - A pointer to the array of column sums of Beta^2.
//rsT - The constraint on the L2 sum of each row of Theta.
//S - Number of samples.
//L - Number of chromosomal locations.
//J - Number of latent features.
/***********************************************************************/
void UpdateTLatL2C(double* pnewT, double* pY, double* pB,
		   double* pBSqCS, double rsT, int S, int L, int J) {

  double *ptldY = new double[S];

  for (int j = 0; j < J; j++) {

    MakeTldY(ptldY,pY,pB,pnewT,j,S,L,J);
    double tldYSqS = 0.0, tldYnrm;
    for (int s = 0; s < S; s++) {
      tldYSqS += ptldY[s]*ptldY[s];
    }
    tldYnrm = sqrt(tldYSqS);

    if (pBSqCS[j]==0) {
      for (int s = 0; s < S; s++) {
	pnewT[j+J*s] = 0.0;
      }
    } else if (tldYnrm <= sqrt(rsT)*(pBSqCS[j])) {
      for (int s = 0; s < S; s++) {
	pnewT[j+J*s] = ptldY[s]/pBSqCS[j];
      }
    } else {
      for (int s = 0; s < S; s++) {
	pnewT[j+J*s] = sqrt(rsT)*ptldY[s]/tldYnrm;
      }
    }

  }
  delete [] ptldY;
}

/************************************************************************/
//The function LatRSS() calculates the residual sum of squares.
//Arguments:
//pY - Pointer to the original data matrix.
//pB - Pointer to the current Beta matrix.
//pT - Pointer to the current Theta matrix.
//S - Number of samples.
//L - Number of chromosomal locations.
//J - Number of latent features.
/************************************************************************/
double LatRSS(double* pY, double* pB, double* pT, int S, int L, int J) {

  double term1 = 0.0;

  for (int lY = 0; lY < L; lY++) {
    for (int sY = 0; sY < S; sY++) {

    double term2 = 0.0;
    for (int j = 0; j < J; j++) {
      term2 += pB[lY+j*L]*pT[j+sY*J];
    }

    term1 += (pY[lY+sY*L] - term2)*(pY[lY+sY*L] - term2);

    }
  }

  return(term1);

}

/************************************************************************/
//The function LatBIC() calculates the BIC.
//Arguments:
//rss - The residual sum of squares.
//pB - Pointer to the current Beta matrix.
//S - Number of samples.
//L - Number of chromosomal locations.
//J - Number of latent features.
/************************************************************************/
double LatBIC(double rss, double* pB, int S, int L, int J) {

  int nparms = 0;
  double bic;

  for (int j = 0; j < J; j++) {
    double *pBcol = new double[L];
    for (int l = 0; l < L; l++) {
      pBcol[l] = pB[l+j*L];
    }
    R_rsort(pBcol,L);
    for (int l = 0; l< L; l++) {
      if ((l==0)&&(pBcol[l]!=0)) {
	  nparms++;
      } else if ((pBcol[l]!=0)&&(pBcol[l]!=pBcol[l-1])) {
	  nparms++;
      }
    }
    delete [] pBcol;
  }

  bic = log(S*L*1.0)*nparms + S*L*1.0*log(rss/(S*L));

  return(bic);

}
