#ifndef FUSED_LASSO_C_SRC_FL_H__

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Print.h>


SEXP L2L1Vit(SEXP obsSeq, SEXP obsWts, SEXP lambda2, 
			 SEXP retPath, SEXP maxSegs, SEXP nSegs, SEXP backPtrs);

SEXP L2L1VitPath(SEXP obsSeq, SEXP lambda2, SEXP retPath, SEXP maxSegs,
				 SEXP segmentVec);

#endif
