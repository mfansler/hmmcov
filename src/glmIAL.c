#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "utility.h"
#include "glm.h"
#include "IAL.h"

/**********************************************************************
 *
 * glmIAL
 *

 
 Input:
 
 family       GLM family (see below)
 link         Link function (see below)
 N            # units
 M            # X variables
 y            y-variable (N-vector)
 X            If M>0, N*M matrix of X variables
 maxit        Maximum number of iterations of IRLS algorithm
 conv         Proportional change in weighted sum of squares residuals to
              declare convergence
 init         If true (non-zero), the iteration starts from initial estimates 
              of fitted values (see below). This option has no effect if
              no iteration is required
 
 Output:
 
 fitted       fitted values 
 resid        working residuals (on linear predictor scale) (N-vector)
 weights      weights (N-vector)
 beta         regression coeficients

 Return
 
 0            convergence
 1            no convergence after maxit iterations
 
 **********************************************************************/

int glmIAL(int* familyR, int *linkR, int *dims, int *nIter,
           double *y, double *offset, double *RX, 
           double *RwX, double *RnX, double *nTotal_binom,  
           double *convR, double *fitted, double *resid,  
           double *weights, double *phi, int *trace, int *n2use,  
           int *w2use, double *b2use, double* b02use, double* Rb0M, 
           double *Rscore, double *score2use, double *delta, 
           double *delta2use, double *tau, double *tau2use, 
           double *wssV2, double *wss0,
           double *prior, double *prop)
{
  double epsilon = 1e-8;       /* Singularity threshold */
  int N, M, maxit, maxitIAL, init, useOffset, n_delta, n_tau;
  int i, j, k, k1, k2, Nu, p_max, n2kpIAL=0, *w2kpIAL;
  int convged = 0, iter = 0, useInitIAL = 0, nReEstimate;
  
  int family = *familyR;
  int link   = *linkR;
  
  if(family > 6 || family < 1){
    Rprintf("family=%d, ", family);
    error("Invalid family!\n");
  }
  
  double conv = *convR;
  double mu, eta, ri, wi, wsum, D, Vmu, wss, wss_last,pi;
  double **score, *z, *wz, xij, **X, **wX, **nX, **b0M;    
  double *b, *v, *Xj2, *Xj, *wXj, BIC=0.0, *resIAL;
  double delta1, tau1, *b2kpIAL, b0=0.0, Xjnorm, *wssV1,penalty=0;
  
  *score2use = DBL_MAX;

  N         = dims[0];
  M         = dims[1];
  maxit     = dims[2];
  maxitIAL  = dims[3];
  init      = dims[4];
  useOffset = dims[5];
  n_delta   = dims[6];
  n_tau     = dims[7];
  p_max     = dims[8];
  nReEstimate = dims[9];
  
  /* Note, each row of X and wX is one covariate */
  reorg(RX,  &X,  M, N);
  reorg(RwX, &wX, M, N);
  k = nReEstimate+1;
  if(k < 1) { k = 1; }
  reorg(RnX, &nX, k, N);  
  reorg(Rscore, &score, n_delta, n_tau);
  reorg(Rb0M,   &b0M,   n_delta, n_tau);

  z  = (double *)Calloc(N, double);  // working y
  wz = (double *)Calloc(N, double);  // weighted working y
  /* variables to be send to function IAL */
  b  = (double *)Calloc(M, double);
  v  = (double *)Calloc(M, double);
  
  wssV1   = (double *)Calloc(maxit, double);
  resIAL  = (double *)Calloc(N, double);
  w2kpIAL = (int *)Calloc(M, int);
  b2kpIAL = (double *)Calloc(M, double);

  /**
   * sum square for each marker, i.e., each column of X
   */
  Xj2 = (double *)Calloc(M, double);
  
  if(*trace > 7){
    Rprintf("y\n");
    Rprint_ve(y, 0, 9);
    Rprintf("X\n");
    Rprint_me(X, 0, 1, 0, 9);
  }
  
  if(*trace>3)
    Rprintf("  glmIAL: family=%d, phi=%f, init=%d\n", family, *phi, init);

  /* ----------------------------------------------------------------
   * lopp across different deltas and taus
   * ----------------------------------------------------------------*/
  
  for(k1=0; k1 < n_delta; k1++){
    
    delta1 = delta[k1];
    
    if(*trace > 5){
      Rprintf("\nk1=%d, delta1=%f\n", k1, delta1);
    }
    
    for(k2=0; k2 < n_tau; k2++){
      
      tau1 = tau[k2];
      
      if(*trace > 5){
        Rprintf("\n  k2=%d, tau1=%f\n", k2, tau1);
      }
      
      /* ----------------------------------------------------------------*
       * by default, initialize mu (fitted) by y itself, with neccesary 
       * modification, e.g., y + 0.1 to avoid log(0) for Poisson family
       * ----------------------------------------------------------------*/
      if (!init) {
        initialize(family, y, fitted, N, nTotal_binom);
      }
      
      /* ----------------------------------------------------------------*
       * Initialize wi (weights) and (standardized) residual 
       * In IRLS, z_i = eta_i + (y_i - mu_i)*(d_eta/d_mu)_i
       *
       * eta_i = linkfun(mu_i)
       * mu_i  = invlink(eta_i)
       *
       * In the following code:
       * ri = (y_i - mu_i)*(d_eta/d_mu)_i
       * wi = Vmu is the weight,  
       *
       * for those invlaid values, set their weight to be 0.
       * ----------------------------------------------------------------*/
      
      wsum = 0.0;
      wss  = 0.0;
      
      for (i=0; i<N; i++) {
    	
	pi=prior[i];
        mu = fitted[i];
        
        if (!muvalid(family, mu)) { 
          wi = ri = 0.0; 
        }else {
          Vmu = varfun(family, mu, *phi);
          
          if (link == family) {
            ri = (y[i] - mu)/Vmu;
            wi = pi*Vmu;
          }else {
            D  = dlink(link, mu);
            ri = D*(y[i] - mu);
            wi = pi/(D*D*Vmu);
          }
          
          wss += wi*ri*ri;
        }
        
        weights[i] = wi;
        resid[i]   = ri;
        if (weights[i] < epsilon) weights[i] = 0.0;
        
        wsum += weights[i];
      }
      
      /* this is the best wss0 we can get if fitted value is initialized by y */
      *wss0 = wss;

      /* ----------------------------------------------------------------*
       * If summation of all weights is too small, stop 
       * ----------------------------------------------------------------*/
      
      if (wsum < epsilon) {
        Rprintf("  glmIAL: summation of all weights are too small!\n");
      }
      
      if(*trace > 5){
        Rprintf("\n  glmIAL: finish initialization, N=%d, M=%d, family=%d\n", 
                N, M, family);
      }
      
      /* ----------------------------------------------------------------*
       * IRLS algorithm  
       * ----------------------------------------------------------------*/
      
      convged  = 0;
      iter     = 0;
      wss_last = 0.0;

      while(iter<maxit && !convged) {
        
        if (*trace > 5) {
          Rprintf("\n  glmIAL: iteration %d: \n", iter);
        }
                
        for (i=0; i<N; i++) {
          /**
           * current estimate of eta + (y-mu)/gradient
           *
           * linkfun(link, fitted[i]) = eta
           * resid[i] = (y-mu)/gradient
           */
          z[i] = linkfun(link, fitted[i]) + resid[i];
        }
        
        if (useOffset) {
          for (i=0; i<N; i++) z[i] -= offset[i];
        }
        
        /**
         * incoporate the weights to z 
         */
        
        for (i=0; i<N; i++) {
          wz[i] = sqrt(weights[i])*z[i];
        }
        
        /**
         * incoporate the weights to X
         * and calculate the scale of X
         */
        
        for (j=0; j<M; j++) {
          Xj     = X[j];
          wXj    = wX[j];
          Xjnorm = 0.0;
          
          for (i=0; i<N; i++) {
            xij     = sqrt(weights[i])*Xj[i];
            wXj[i]  = xij;
            Xjnorm += xij*xij;
          }
          
          Xj2[j] = Xjnorm;
        }
        
        if (iter==0) {
          useInitIAL=0;
        }else {
          useInitIAL=1;
        }

        IAL(wz, wX, nX, &b0, b, v, delta1, tau1, dims, Xj2, conv, resIAL,  
            &BIC, &n2kpIAL, w2kpIAL, b2kpIAL, useInitIAL, *trace, nReEstimate, weights,prop);
        
//        if (n2kpIAL == 0) { break; }
        
        /**
         * update the fitted values and weights  
         */
        
        wss = 0.0;
        Nu  = 0;
        
        for (i=0; i<N; i++) {
          
          eta = b0;
          for (k=0; k<n2kpIAL; k++) {
            j    = w2kpIAL[k];
            eta += X[j][i]*b2kpIAL[k];
          }
           
          if (useOffset) { eta += offset[i]; }
          
          pi=prior[i];
          mu = invlink(link, eta);
          
          fitted[i] = mu;
          
          if (weights[i]<=0.0) {
            wi = ri = 0.0;
          } else {
            
            Vmu = varfun(family, mu, *phi);
            Nu ++;
            
            if (link == family) {
              ri = (y[i] - mu)/Vmu;
              wi = pi*Vmu;
            }else {
              D  = dlink(link, mu);
              ri = D*(y[i] - mu);
              wi = pi/(D*D*Vmu);
            }
            if (wi < epsilon) wi = 0.0;

            wss += wi*ri*ri;
            
          }
          
          weights[i] = wi;
          resid[i]   = ri;
        }
        
        if(*trace > 9){
          if (n2kpIAL > 0) {
            Rprintf("w2kpIAL: ");
            Rprint_vi(w2kpIAL, 0, n2kpIAL-1);
            Rprintf("b2kpIAL: ");
            Rprint_ve(b2kpIAL, 0, n2kpIAL-1);        
          }
        }
        
        if(wss > 1/epsilon){
          if(*trace > 5){
            Rprintf("  glmIAL: huge wss (%,3e)!\n", wss);
          }
          
//          break;
        }
        
        if (*trace > 5) {
          Rprintf("    n2kpIAL=%d, wss=%.3e\n", n2kpIAL, wss);
        }
        
        wssV1[iter] = wss;
        convged  = (Nu<=0) || (fabs(wss-wss_last)/(wss_last + 0.1) < conv);
        wss_last = wss;
        iter++;
                
      }
      
      if(*trace > 5){
        Rprintf("    n2kpIAL=%d, BIC=%.3e\n", n2kpIAL, BIC);
      }
      
      if(convged==0 /*|| n2kpIAL > p_max*/ & *trace>5){	
	Rprintf("glmIAL: Failure to Converge\n");
      }

      //Hopefull this will not result in an error, but want to report even if just intercept
      if(n2kpIAL==0){
//        continue;
      }
      
      if(*trace > 5){
        if (n2kpIAL > 0) {
          Rprintf("w2kpIAL: ");
          Rprint_vi(w2kpIAL, 0, n2kpIAL-1);
          Rprintf("b2kpIAL: ");
          Rprint_ve(b2kpIAL, 0, n2kpIAL-1);          
        }
      }
      
      score[k1][k2] = BIC;
      b0M[k1][k2]   = b0;

      if (BIC < *score2use) {
        *score2use = BIC;
        *delta2use = delta1;
        *tau2use   = tau1;
        *n2use     = n2kpIAL;
        *b02use    = b0;
        *nIter     = iter;

        for (i=0; i < n2kpIAL; i++) {
          w2use[i] = w2kpIAL[i];
          b2use[i] = b2kpIAL[i];
	  			//penalty += *prop*(1 + *delta2use)*log(fabs(b2use[i]) + *tau2use);
	  			//*score2use = penalty;
        }
        
        for (i=n2kpIAL; i < p_max; i++) {
          w2use[i] = -9;
          b2use[i] = 0.0;
        }
        
        for (i=0 ; i < iter; i++) {
          wssV2[i] = wssV1[i];
        }
        
        if(iter < maxit){
          for (i=iter ; i < maxit; i++) {
            wssV2[i] = 0.0;
          }
        }
        
      }

    }
  }
  
  if(*trace > 5){
    Rprintf("\nn2use = %d, score2use=%.2e\n", *n2use, *score2use);
    Rprintf("score: \n");
    Rprint_me(score, 0, n_delta-1, 0, n_tau-1);

    if (*n2use > 0) {
      Rprintf("w2use: \n");
      Rprint_vi(w2use, 0, *n2use-1);
      Rprintf("b2use: \n");
      Rprint_ve(b2use, 0, *n2use-1);          
    }
  }

  Free(z);
  Free(wz);
  Free(Xj2);
  Free(b);
  Free(v);
  Free(wssV1);
  Free(w2kpIAL);
  Free(b2kpIAL);
  Free(resIAL);

  return(*n2use);
}
