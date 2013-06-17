/*
 *  IAL.c
 *
 *  Created by Wei Sun on Sep 8 2010.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "IAL.h"
#include <R_ext/Applic.h>

/*********************************************************************
 *
 * IAL
 *
 * The Iterative Adaptive Lasso
 *
 
   y:       response varaible
   X:       covaraite marix, one row for one covariate
   b0:      intercept
   b:       regression coefficient except intercept
   v:       a set of (combination of) penalization parameters
   delta:   tuning parameters delta
   tau:     tuning parameters tau
   dims:    dimensions
   conv:    a small positive number of to call convergence
   resid:   working residuals, has been initialized by y
   BIC:     BIC score
   Xj2:     sum square of Xj
   n2kp:    number of covaraites seleted by IAL
 
 *********************************************************************/

void IALtest(double* y, double** X, double** Xnew, double* b0,   
         double* b, double* v, double delta1, double tau1,    
         int* dims, double *Xj2, double conv, double* resid, 
         double* BIC, int* n2kp, int* w2kp, double* b2kp, 
         int initIAL, int trace, int nReEstimate, 
	 double *weights, double *prop)
{
  int i, j, k, w, n, p, pp, maxitIAL, converged;
  double bj_bar, v0, v0_last, sj, bT, delta2, ave, bj, *Xj;
  
  n = dims[0];
  p = dims[1];
  maxitIAL = dims[3];
    
  /* variables used for the least squares at the end */
  double *nX_j, *X_jk, yi;
  
  /* variables used for the function dqrls 
   *
   F77_CALL(dqrls)(Xnew, n, k, y, ny, tol, coef, resid,
   qty, rank, jpvt, qraux, work);
   **/
  
  int ny=1, rank=1, *jpvt;
  double tol=1e-7, *qty, *coef, *qraux, *work;
  
  qty   = (double *)Calloc(n, double);
  qraux = (double *)Calloc(nReEstimate+1, double);
  coef  = (double *)Calloc(nReEstimate+1, double);
  work  = (double *)Calloc(2*(nReEstimate+1), double);
  jpvt  = (int *)Calloc(nReEstimate+1, int);

  if(trace > 5){
    Rprintf("  IAL: n=%d, p=%d, maxitIAL=%d\n", n, p, maxitIAL);
  }
    
  if(trace > 11){
    Rprintf("y\n");
    Rprint_ve(y, 0, 9);
    Rprintf("resid\n");
    Rprint_ve(resid, 0, 9);    
    Rprintf("X\n");
    Rprint_me(X, 0, 2, 0, 9);
  }
  
  /* IAL calculation for given delta and tau */
  delta2 = *prop*(delta1 + 1.0);
  
  /**
   * step 1. Initialization
   */
  
  if (!initIAL) {
    for(j=0; j<p; j++){ b[j] = 0.0; }
  }
  
  *b0 = 0.0;
  
  /**
   * step 2 Update v0, b0, and resid
   */
  
  v0_last = 0.0;
  
  ave = 0.0;
  v0  = 0.0;
  
  for(i=0; i<n; i++){ 
    yi   = y[i];
    ave += yi;
    v0  += yi*yi;
    resid[i] = yi;
  }
  ave = ave/n;
  v0  = v0/n - ave*ave;
  
  if (fabs(ave) > 1e-16) {
    for(i=0; i<n; i++){ resid[i] -= ave; }
    *b0 = *b0 + ave;
  }
  
  /**
   * step 3 Update v[j], i.e. kappa_j
   */
  for(j=0; j<p; j++){
    v[j] = (fabs(b[j]) + tau1)/delta2;
  }          
    
  for(w=0; w<maxitIAL; w++){
    
    /**
     * step 4 Update b[j]
     */
    for(j=0; j<p; j++){
      Xj = X[j];
      bj = b[j];
      
      /* remove the effect of Xj from the residual*/
      
      bj_bar = 0.0;

      if (fabs(bj) > 1e-16) {
        for(i=0; i<n; i++){ 
          resid[i] += Xj[i]*bj; 
          bj_bar += Xj[i]*resid[i];
        }            
      }else {
        for(i=0; i<n; i++){
          bj_bar += Xj[i]*resid[i];
        }
      }

      bj_bar /= Xj2[j];
      sj      = v0/Xj2[j];
      
      bT = sj/v[j];
      
      if(bj_bar > bT){
        bj = bj_bar - bT;
      }else if(bj_bar < -bT){
        bj = bj_bar + bT;
      }else{
        bj = 0.0;
      }

      /* add the effect of Xj back into the residual */
      if (fabs(bj) > 1e-16) {
        for(i=0; i<n; i++){ resid[i] -= Xj[i]*bj; }
      }
      
      b[j] = bj;
    }

    /**
     * step 1. Update b0 and resid
     * actually this step may not be neccesary
     * this step needs to be looked at again.  Now, resid is y-sum(xb_new)- b0_old 
     */
    
    ave = mean(resid, 0, n-1);
    
    if(fabs(ave) > 1e-16){
      for(i=0; i<n; i++){ resid[i] -= ave; }
      
      *b0 = *b0 + ave;
    }

    /*calculate ybar-sum(bjxbarj)
    *b0 = mean(resid, 0, n-1);
    for(j=0; j<p; j++){
      Xj = X[j];
     bj = b[j];
      *b0=*b0 - bj * mean(Xj, 0, n-1);
    }

    for(i=0; i<n; i++){ resid[i] -= *b0; }*/

    
    /**
     * step 2 Update v0, i.e. sigma^2
     */
    
    v0 = 0.0;
    for(i=0; i<n; i++){ 
      yi  = resid[i];
      v0 += yi*yi; 
    }
    v0 /= n;
    
    /**
     * step 3 Update v[j], i.e. kappa_j
     */
    for(j=0; j<p; j++){
      v[j] = (fabs(b[j]) + tau1)/delta2;
    }          
        
    /**
     * check convergence
     */          
    converged = (fabs(v0 - v0_last)/(v0_last + 0.1) < conv);

    if(trace > 9){
      Rprintf("      w=%d, v0=%f, converged=%d\n", w, v0, converged);
    }
    
    if(converged){ break; }
    v0_last = v0;
  }

  k = 0;
  
  for(j=0; j < p; j++){
    if(fabs(b[j]) > 1e-10){ 
      w2kp[k] = j;
      b2kp[k] = b[j];
      k++; 
    }
  }
  /*
    Rprintf("updated v0=%f\n", v0);
    Rprintf("b0=%f\n", *b0);
    Rprintf("b2kp=: ");
    Rprint_ve(b, 0, p-1);
 */   

  /**
   * Use ordinary least squares refit the model and then estimate BIC
   * if the number of covariates selected is relatively small
   */
  
  if(k > 0 && k < nReEstimate){
    
    if (trace > 7) {
      Rprintf("reestimate coefficients\n");
      Rprintf("k=%d, nReEstimate=%d\n", k, nReEstimate);
    }

    /* Now the dimension of Xnew is (k+1)*n */
    
    nX_j = Xnew[0];
    for (i=0; i<n; i++) {
      nX_j[i] = sqrt(weights[i]);
    }
    jpvt[0] = -1;
    
    for (j=0; j<k; j++) {
      jpvt[j+1] = j;
      
      nX_j = Xnew[j+1];
      X_jk = X[w2kp[j]];
      
      for (i=0; i<n; i++) { nX_j[i] = X_jk[i]; }
    }
   
    ny   = 1;
    pp   = k+1;
    rank = k+1;
    
    F77_CALL(dqrls)(Xnew[0], &n, &pp, y, &ny, &tol, coef, resid,
                    qty, &rank, jpvt, qraux, work);
    
    if (trace > 7) {
      Rprintf("result of F77_CALL(dqrls): rank=%d, jpvt\n", rank);
      Rprint_vi(jpvt, 0, k);
      
      Rprintf("coef\n");
      Rprint_ve(coef, 0, k);
    }

    for (j=0; j<=k; j++) { 
      if (jpvt[j] < 0) {
        *b0 = coef[j];
      }else {
        b2kp[jpvt[j]] = coef[j]; 
      }
    }
    
    v0 = 0.0;
    for(i=0; i<n; i++){ v0 += resid[i]*resid[i]; }
    v0 /= n;
    
   
  }
  
  *BIC  = n*log(v0) + k*log(n);
  *n2kp = k;
  
  Free(qty);
  Free(qraux);
  Free(coef);
  Free(work);
  Free(jpvt);

}
