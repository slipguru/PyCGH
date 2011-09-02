/*************************************************************************/
//This file contains the C function for running the Latent Feature Model
//under the L2 constraint.
/*************************************************************************/

#include <cmath>
#include <R.h>
#include <Rinternals.h>
#include "gen_lat_func.h"

extern "C" {


  SEXP LatL2C(SEXP Y, SEXP nF, SEXP inB, SEXP inT, SEXP lam1, SEXP lam2,
	      SEXP thresh, SEXP maxiter, SEXP maxiterB, SEXP maxiterT,
	      SEXP sT) {

    int nProt = 0;
    double *pY = REAL(Y), *pinB = REAL(inB), *pinT = REAL(inT),
      rlam1 = REAL(lam1)[0], rlam2 = REAL(lam2)[0],
      rthresh = REAL(thresh)[0], rsT = REAL(sT)[0];
    int imaxiter = INTEGER(maxiter)[0], imaxiterB = INTEGER(maxiterB)[0],
      imaxiterT = INTEGER(maxiterT)[0];
    R_len_t S = ncols(Y), L = nrows(Y), J = INTEGER(nF)[0];

    //Initializing Beta.
    SEXP newB = R_NilValue;
    PROTECT(newB = allocMatrix(REALSXP,L,J));
    nProt++;
    double *pnewB = REAL(newB), *poldB = new double[L*J];
    CopyAtoB(pinB,poldB,L*J);
    CopyAtoB(pinB,pnewB,L*J);
    //Initializing Theta.
    SEXP newT = R_NilValue;
    PROTECT(newT = allocMatrix(REALSXP,J,S));
    nProt++;
    double *pnewT = REAL(newT), *poldT = new double[J*S];
    CopyAtoB(pinT,poldT,J*S);
    CopyAtoB(pinT,pnewT,J*S);

    double *perrBs = new double[imaxiter+1],
      *perrTs = new double[imaxiter+1], BSqS = SqTotSum(poldB,L*J);
    perrBs[0] = rthresh +1.0;
    perrTs[0] = rthresh + 1.0;
    int niter = 0;

    while (((perrBs[niter]>rthresh)||(perrTs[niter]>rthresh))&&
	   (niter<imaxiter)&&(BSqS!=0.0)) {
      TLatL2C(pnewT,pY,pnewB,rthresh,imaxiterT,rsT,S,L,J);
      perrTs[niter+1] = MatErr(pnewT,poldT,J*S,rthresh);
      CopyAtoB(pnewT,poldT,J*S);
      BC(pnewB,pY,pnewT,rlam1,rlam2,rthresh,imaxiterB,S,L,J);
      perrBs[niter+1] = MatErr(pnewB,poldB,L*J,rthresh);
      CopyAtoB(pnewB,poldB,L*J);
      BSqS = SqTotSum(poldB,L*J);
      niter++;
    }

    SEXP Rniter = R_NilValue, rss = R_NilValue, bic = R_NilValue;
    PROTECT(rss = allocVector(REALSXP,1));
    nProt++;
    REAL(rss)[0] = LatRSS(pY,pnewB,pnewT,S,L,J);
    PROTECT(bic = allocVector(REALSXP,1));
    nProt++;
    REAL(bic)[0] = LatBIC(REAL(rss)[0],pnewB,S,L,J);
    PROTECT(Rniter = allocVector(INTSXP,1));
    nProt++;
    INTEGER(Rniter)[0] = niter;

    //The results.
    SEXP res = R_NilValue, resNames = R_NilValue, resClass = R_NilValue;
    //The list of output variables.
    PROTECT(res = allocVector(VECSXP,5));
    nProt++;
    SET_VECTOR_ELT(res,0,newB);
    SET_VECTOR_ELT(res,1,newT);
    SET_VECTOR_ELT(res,2,Rniter);
    SET_VECTOR_ELT(res,3,rss);
    SET_VECTOR_ELT(res,4,bic);
    //The list of output variable names.
    PROTECT(resNames =allocVector(STRSXP,5));
    nProt++;
    SET_STRING_ELT(resNames,0,mkChar("Beta"));
    SET_STRING_ELT(resNames,1,mkChar("Theta"));
    SET_STRING_ELT(resNames,2,mkChar("niter"));
    SET_STRING_ELT(resNames,3,mkChar("rss"));
    SET_STRING_ELT(resNames,4,mkChar("bic"));
    //Setting the names to the output variables.
    setAttrib(res,R_NamesSymbol,resNames);
    //Setting the class to the output list.
    PROTECT(resClass = allocVector(STRSXP,1));
    nProt++;
    SET_STRING_ELT(resClass,0,mkChar("FLLat"));
    classgets(res,resClass);

    delete [] perrBs;
    delete [] perrTs;
    delete [] poldB;
    delete [] poldT;
    UNPROTECT(nProt);
    return(res);

  }


}
