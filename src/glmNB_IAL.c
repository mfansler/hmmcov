
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
 * glmIAL_one
 * comparing with glmIAL, glmIAL_one only take one combination of
 * penalization parameter
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

void glmIAL_one(int family, int link, int *dims, double *y, 
                double *offset, double **X, int *succeed,
                double **wX, double **nX, double *nTotal_binom,  
                double conv, double *fitted, double *resid,  
                double *weights, double *phi, int *trace, 
                int *n2kpIAL, double *b0, double *BIC, 
                double delta1, double tau1, double *z, 
                double *wz, double *b, double *v, double *resIAL, 
                int *w2kpIAL, double *b2kpIAL, double *Xj2,
		double *prior, double *prop)

{
  double epsilon = 1e-8;       /* Singularity threshold */
  int N, M, maxit, maxitIAL, init, useOffset;
  int i, j, k, k1, k2, Nu, p_max, n_delta, n_tau;
  int convged = 0, iter = 0, useInitIAL = 0, nReEstimate;
  
  if(family > 6 || family < 1){
    Rprintf("family=%d, ", family);
    error("Invalid family!\n");
  }
  
  double mu, eta, ri, wi, wsum, D, Vmu, wss, wss_last, pi;
  double xij, *Xj, *wXj, Xjnorm;
  
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
  
  
  if(*trace > 7){
    Rprintf("y\n");
    Rprint_ve(y, 0, 9);
    Rprintf("X\n");
    Rprint_me(X, 0, 1, 0, 9);
  }
  
  if(*trace>3)
    Rprintf("  glmIAL: family=%d, phi=%f, init=%d\n", family, *phi, init);
  
  useInitIAL=0;
  
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
    }
    
    weights[i] = wi;
    resid[i]   = ri;
    if (weights[i] < epsilon) weights[i] = 0.;
    
    wsum += weights[i];
  }
  
  
  /* ----------------------------------------------------------------*
   * If summation of all weights is too small, stop 
   * ----------------------------------------------------------------*/
  
  if (wsum < epsilon) {
    if(*trace>2)
      Rprintf("  glmIAL_one: summation of all weights are too small!\n");
    
    //maxit    = -1;
    //*succeed = 0;
  }else if(*trace > 5){
    Rprintf("\n  glmIAL_one: finish initialization, N=%d, M=%d, family=%d\n", 
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
      wz[i]     = sqrt(weights[i])*z[i];
      resIAL[i] = wz[i];
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
    
		if(i > 0) useInitIAL = 1;

    IAL(wz, wX, nX, b0, b, v, delta1, tau1, dims, Xj2, conv, resIAL,  
        BIC, n2kpIAL, w2kpIAL, b2kpIAL, useInitIAL, *trace, nReEstimate, weights, prop);
        
    /**
     * update the fitted values and weights  
     */
    
    wss = 0.0;
    Nu  = 0;
    
    for (i=0; i<N; i++) {
      
      eta = *b0;
      for (k=0; k < *n2kpIAL; k++) {
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
    
    if(wss > 1/(epsilon*epsilon*epsilon)){
      if(*trace > 5){
        Rprintf("  glmIAL: huge wss (%,3e)!\n", wss);
      }
      
      convged  = 0;
      *succeed = 0;
      //break;
    }
    
    if (*trace > 5) {
      Rprintf("    n2kpIAL=%d, wss=%.3e\n", *n2kpIAL, wss);
    }
    
    convged  = (Nu<=0) || (fabs(wss-wss_last)/(wss_last + 0.1) < conv);
    wss_last = wss;
    iter++;
    
  }
  
  if(*trace > 3){
    Rprintf("    end of glmIAL_one: n2kpIAL=%d, BIC=%.3e\n", *n2kpIAL, *BIC);
  }
}

/**********************************************************************
 *
 * main function of glmNB
 *
 **********************************************************************/

void glmNB_IAL(int *dims, double *y, double *offset,double *RX, double *convR,  
							 double *scoreTestP, int *trace, double *delta, double *tau,int *nIter, 
							 double *phi, int *n2use, int *w2use,  double *b2use, double *b02use, 
							 double *Rb0M, double *Rscore, double *score2use,  double *delta2use, double *tau2use, 								 double *likelihood, int *family, double *prior, double *prop, double *fitted)
{
  int N, M, maxit, maxitIAL, init, useOffset, succeed=1;
  int i, j, k, k1, k2, kit, iter = 0, n_delta, n_tau, p_max;
  int fam0=0, linkR=0, n2kpIAL, *w2kpIAL, nReEstimate;

  double conv = *convR, eta=0.0, b0=0.0, BIC=0.0;
  double nTotal_binom=0.0;  /* for binomial link, NOT useful here */
  double del=1.0, Lm=0.0, Lm0=0.0, phi0=0.0, delta1, tau1;
  double Dscore, scoreNum, scoreDen, scorePval, yi, mui;
  
  double *z, *wz, **X, **wX, **nX;
  double *b, *v, *Xj2, *resIAL, *b2kpIAL, *lkhood;
  double /**fitted,*/ *resid, *weights, penalty; 
  
  /* convergence indicator for phi_ml 
   * if cvPhi = 0, NB model is OK. 
   * if cvPhi = 1, suggest we need to use Poisson
   * if cvPhi = 2, suggest we need to use ZINB
   */
  int cvPhi;

  *n2use = 0;
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
  
  X    = (double **)Calloc(M, double*);
  X[0] = RX;
  for (j=1; j<M; j++) { X[j] = RX + j*N; }
  
  /* wX means weighted X */
  wX    = (double **)Calloc(M, double*);
  wX[0] = (double *) Calloc(M*N, double);
  for(j=1; j<M; j++){ wX[j] = wX[0] + j*N; }

  nX    = (double **)Calloc(nReEstimate+1, double*);
  nX[0] = (double *) Calloc((nReEstimate+1)*N, double);
  for(j=1; j<=nReEstimate; j++){ nX[j] = nX[0] + j*N; }
  
  //fitted  = (double *)Calloc(N, double);
  resid   = (double *)Calloc(N, double);
  weights = (double *)Calloc(N, double);

  z  = (double *)Calloc(N, double);  // working y
  wz = (double *)Calloc(N, double);  // weighted working y
  
  /* variables to be send to function IAL */
  b  = (double *)Calloc(M, double);
  v  = (double *)Calloc(M, double);
  
  resIAL  = (double *)Calloc(N, double);
  w2kpIAL = (int *)Calloc(M, int);
  b2kpIAL = (double *)Calloc(M, double);
  lkhood  = (double *)Calloc(maxit, double);  // working likelihood

  /* sum square for each marker, i.e., each column of X */
  Xj2 = (double *)Calloc(M, double);
  
  if(*trace > 3) 
    Rprintf("\n  glmNB: N=%d, M=%d, maxit=%d, init=%d, useOffset=%d\n", 
            N, M, maxit, init, useOffset);
  
  /* ----------------------------------------------------------------
   * lopp across different deltas and taus
   * ----------------------------------------------------------------*/
  
  for(k1=0; k1 < n_delta; k1++){
    
    delta1 = delta[k1];
    
    if(*trace > 3){
      Rprintf("\nk1=%d, delta1=%f\n", k1, delta1);
    }
    
    for(k2=0; k2 < n_tau; k2++){

      succeed = 1;
      tau1 = tau[k2];
      
      iter = 0;

      if(*trace > 3){
        Rprintf("\n  k2=%d, tau1=%f\n", k2, tau1);
      }
      
      if(*phi < 0){
	
		    /* Initial fit */
	      fam0  = POISSON;
	      linkR = LOG;
	      dims[4] = 0; /* do not use initial values */
	
	      glmIAL_one(fam0, linkR, dims, y, offset, X, &succeed,
	                 wX, nX, &nTotal_binom, conv, fitted,  
	                 resid, weights, phi, trace, &n2kpIAL, &b0, 
	                 &BIC, delta1, tau1, z, wz, b, v, resIAL, 
	                 w2kpIAL, b2kpIAL, Xj2, prior, prop);

	      if(!succeed & *trace>5){ Rprintf("glmNBIAL_one: Failure to converge");}
  		   
      
			/*
       WELL,  we may alow the case that n2kpIAL = 0 for isoform mapping
       since this indicate that there is only one isoform.
       
      if (n2kpIAL==0) {
        if(*trace > 3){
          Rprintf("\n  glmNB_IAL: no covariate is selected by baseline model\n");
        }
        continue;
      }
      */
      
      /* update fitted value *

      	for (i=0; i<N; i++) {
        
        	eta = b0;
        	for (k=0; k < n2kpIAL; k++) {
        	  j    = w2kpIAL[k];
        	  eta += X[j][i]*b2kpIAL[k];
        	}
        
        	if (useOffset) { eta += offset[i]; }
        
       	 fitted[i] = invlink(linkR, eta);
      	}
  	    
  	    /* test for overdispersion by Dean's Score test *
  	    scoreNum = 0.0;
  	    scoreDen = 0.0;
      
	      for (i=0; i<N; i++) {
        	yi  = y[i];
        	mui = fitted[i];
        	scoreNum += (yi - mui)*(yi - mui) - yi;
        	scoreDen += mui*mui;
      	}
      
      	Dscore = scoreNum/sqrt(2.0*scoreDen);
      
      	/**
      	 * double pnorm(double x, double mu, double sigma, int lower_tail, int give_log);
      	 *
      	scorePval = pnorm(Dscore, 0.0, 1.0, 0, 0);
      
      	if(*trace > 3) 
      	  Rprintf("\n  overdispersion Dean's score = %.2e, p-value = %.2e\n\n", 
                Dscore, scorePval);
      
     	 /* use Poisson model *
      		if(scorePval > *scoreTestP){
        
        	//if(n2kpIAL > p_max){ continue; }

        	kit         = k1*n_tau + k2;
        	Rscore[kit] = BIC;
        	Rb0M[kit]   = b0;
        
        	if (BIC < *score2use) {
        	  *score2use = BIC;
        	  *delta2use = delta1;
        	  *tau2use   = tau1;
        	  *n2use     = n2kpIAL;
        	  *b02use    = b0;
        	  *nIter     = iter;
        	  *family    = fam0;
          
        	  for (i=0; i < n2kpIAL; i++) {
        	    w2use[i] = w2kpIAL[i];
        	    b2use[i] = b2kpIAL[i];
        	  }
        	  
        	  for (i=n2kpIAL; i < p_max; i++) {
        	    w2use[i] = -9;
        	    b2use[i] = 0.0;
        	  }
        	}
        
        	continue;
      	}
        
      	/* use Negative binomial model *
      	fam0  = NB;
      
      	/**
       	* calculate phi by MLE, without initial values of phi
       	*
      
      	cvPhi = phi_ml(y, fitted, N, maxit, conv, phi, 0, *trace, prior);
      
      	if(cvPhi==0){
        
        	if(*trace > 3) 
        	  Rprintf("\n  initial estimate of phi: %e\n", *phi);
        
      	}else if (cvPhi==1){
        
        	if(*trace > 3) 
        	  Rprintf("\n  Choose Poisson family due to small phi:%e\n", *phi);
        
       	/* if(n2kpIAL > p_max){ continue; }
        
        	kit         = k1*n_tau + k2;
        	Rscore[kit] = BIC;
        	Rb0M[kit]   = b0;
        
        	if (BIC < *score2use) {
        	  *score2use = BIC;
        	  *delta2use = delta1;
        	  *tau2use   = tau1;
        	  *n2use     = n2kpIAL;
        	  *b02use    = b0;
        	  *nIter     = iter;
        	  *family    = fam0;
          
        	  for (i=0; i < n2kpIAL; i++) {
        	    w2use[i] = w2kpIAL[i];
        	    b2use[i] = b2kpIAL[i];
        	  }
          
        	  for (i=n2kpIAL; i < p_max; i++) {
        	    w2use[i] = -9;
        	    b2use[i] = 0.0;
        	  }
        	}
        
        	continue;*
        
      	}else if(cvPhi==2){
        
        	if(*trace > 3) 
        	  Rprintf("\n  The overdispersion parameter is too large: phi=%e\n", *phi);
        
        	//continue;
        
      	}else { /* estimation of phi fail to converge *
        	if(*trace > 3)
        	  Rprintf("\n  glmNB_IAL: fail to converge in phi_ml: phi=%e\n", *phi);
        
        	//continue;
      	}
			}
      /**
       * iterative updates
       */

			cvPhi = phi_ml(y, fitted, N, maxit, conv, phi, 0, *trace, prior);
     	/*fam0  = NB;
	    linkR = LOG;
      del  = 1.0;
      Lm   = loglik_NB(N, *phi, fitted, y, prior);
			Lm0  = Lm + 1.0;
			iter = 0;*/
		}//else{
     	fam0  = NB;
	    linkR = LOG;
      Lm   = loglik_NB(N, *phi, fitted, y, prior);
			Lm0  = Lm;
			del = 0;
      iter = 0;
		//}


			dims[4] = 1; 
			fam0  = NB;
	    linkR = LOG;
      glmIAL_one(fam0, linkR, dims, y, offset, X, &succeed, 
                   wX, nX, &nTotal_binom, conv, fitted,  
                   resid, weights, phi, trace, &n2kpIAL, &b0, 
                   &BIC, delta1, tau1, z, wz, b, v, resIAL, 
                   w2kpIAL, b2kpIAL, Xj2, prior, prop);

		
			if(*trace>6) Rprintf("finish first and only glm iteration\n");

      
      while (iter < maxit && fabs(Lm0 - Lm) + fabs(del) > conv) {
        //Rprintf("%f\n" ,fam0);
        if (*trace > 3) {
          Rprintf("\n  iteration %d in glmNB_IAL\n", iter);
        }
        
        dims[4] = 1; /* use initial values */
        
        glmIAL_one(fam0, linkR, dims, y, offset, X, &succeed, 
                   wX, nX, &nTotal_binom, conv, fitted,  
                   resid, weights, phi, trace, &n2kpIAL, &b0, 
                   &BIC, delta1, tau1, z, wz, b, v, resIAL, 
                   w2kpIAL, b2kpIAL, Xj2, prior, prop);

      
        if(!succeed & *trace>5){ Rprintf("glmNBIAL_one: Failure to converge");}

        /*
         WELL,  we may alow the case that n2kpIAL = 0 for isoform mapping
         since this indicate that there is only one isoform.
         
         if (n2kpIAL==0) {
           if(*trace > 3){
             Rprintf("\n  glmNB_IAL: find 0 covariate in glmIAL_one\n");
           }
           break;
         }
         */
        
        if(*trace > 3){
          Rprintf("\n  n2kpIAL=%d, b0=%f, BIC=%.2e\n", n2kpIAL, b0, BIC);
          
          if (n2kpIAL > 0) {
            Rprintf("  w2kpIAL: ");
            Rprint_vi(w2kpIAL, 0, n2kpIAL-1);
            Rprintf("  b2kpIAL: ");
            Rprint_ve(b2kpIAL, 0, n2kpIAL-1);          
          }
        }
        
        /* update fitted value */
        
        for (i=0; i<N; i++) {
          
          eta = b0;
          for (k=0; k < n2kpIAL; k++) {
            j    = w2kpIAL[k];
            eta += X[j][i]*b2kpIAL[k];
          }
          
          if (useOffset) { eta += offset[i]; }
          
          fitted[i] = invlink(linkR, eta);
        }
        
        phi0  = *phi;
        
        cvPhi = phi_ml(y, fitted, N, maxit, conv, phi, 1, *trace, prior);
        
        if(cvPhi==0){
          if(*trace > 3) 
            Rprintf("\n  finish phi_ml, cvPhi=%d, phi=%e\n", cvPhi, *phi);
          
        }else if(cvPhi==1){
          if(*trace > 1) 
            Rprintf("  glmNB_IAL: overdispersion is too small: phi=%e\n", *phi);
          
          //succeed = 0;
          //break;
        }else if(cvPhi==2){
          if(*trace > 1) 
            Rprintf("  glmNB_IAL: too large overdispersion: phi=%e\n", *phi);
          
          //succeed = 0;
          //break;
        }else {
          if(*trace > 1) 
            Rprintf("  glmNB_IAL: fail to converge in phi_ml\n");
          
          //succeed = 0;
          //break;
        }
        
        del = phi0 - *phi;
        Lm0 = Lm;
        Lm  = loglik_NB(N, *phi, fitted, y, prior);
        
        if (*trace > 3) {
          Rprintf("\n  Phi(%d) = %.2e, Lm=%f, Lm0=%f, del=%f\n\n", 
                  iter, *phi, Lm, Lm0, del);
        }
        lkhood[iter] = Lm;
        iter++;
      }
      
      //if (!succeed) { continue; }
      
      if(iter == maxit) {
        if (*trace > 1) {
          Rprintf("\n  delta=%f, tau=%f\n", delta1, tau1);
          Rprintf("  glmNB_IAL: Alternation limit reached: iter=%d\n", iter);
        }
      }

      //if(n2kpIAL > p_max){ continue; }
      
      kit         = k1*n_tau + k2;
      Rscore[kit] = BIC;
      Rb0M[kit]   = b0;

      if (BIC < *score2use) {
        *score2use = BIC;
        *delta2use = delta1;
        *tau2use   = tau1;
        *n2use     = n2kpIAL;
        *b02use    = b0;
        *nIter     = iter;
        *family    = fam0;
        
        penalty = 0.0;
        
        for (i=0; i < n2kpIAL; i++) {
          w2use[i] = w2kpIAL[i];
          b2use[i] = b2kpIAL[i];
          penalty += *prop*(1 + *delta2use)*log(fabs(b2use[i]) + *tau2use);
        }
        
        for (i=n2kpIAL; i < p_max; i++) {
          w2use[i] = -9;
          b2use[i] = 0.0;
        }

        for (i=0; i < iter; i++) {
          likelihood[i] = lkhood[i];
        }
                
        for (i=iter; i < maxit; i++) {
          likelihood[i] = 0.0;
        }
        
      }
      
    }// end of loop for tau
  }// end of loop for delta
   
  
  Free(z);
  Free(wz);
  Free(b);
  Free(v);
  Free(resIAL);
  Free(w2kpIAL);
  Free(b2kpIAL);
  Free(Xj2);
  Free(lkhood);
  
  Free(X);
  Free(X[0]); //new
  Free(wX[0]);
  Free(wX);
  Free(nX[0]);
  Free(nX);
  
  //Free(fitted);
  Free(resid);
  Free(weights);
}

