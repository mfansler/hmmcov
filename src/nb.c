#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "utility.h"
#include "lbfgsb1.h"

/*double digamma(double *x)
{
  const double h=1.0e-6;
  return ((lgamma(*x + h) - lgamma(*x - h))/(2 * h));
}
*/

void NegbinGradtest(double* p, double* y, double *w, double* ee, double* x, int* nrowx, int* ncolx, double* grad)
{
// first derivative: negative log-likelihood of censored zero-inflated negative binomial
	int i,j, n=*nrowx, colx = *ncolx;
	double loglambda,lambda,logtau, tau_yt,tau;
	for(j=0; j<colx;j++)  grad[j] = 0; // initialize grad
	tau =p[colx];
	
	for(i=0; i<n; i++)
	{
	 loglambda =0;
	 for(j=0; j<colx;j++) loglambda += x[i + j*n]*p[j];
	 lambda = exp(loglambda)*ee[i];
	 loglambda = log(lambda);
	 for(j=0; j<colx;j++)
	    	 //grad[j]         += x[i + j*n]*(1-(tau + y[i])/(lambda + tau))*tau;
          grad[j]         += w[i]*(x[i + j*n]*tau*(lambda-y[i])/(lambda + tau));
          tau_yt = tau + y[i];
          grad[colx] +=   w[i]*(-1  + (tau_yt)/(lambda + tau) + log((lambda + tau)/tau)
																+ digamma(tau) - digamma(tau_yt));
                               //+ digamma(&tau) - digamma(&tau_yt));
   }
  grad[colx]  = tau*grad[colx];
}

void NegbinNLLtest(double* p, double* y, double *w, double* ee, double* x, int* nrowx, int* ncolx, double* nll)
{
// negative log-likelihood of negative binomial
	int i, j,  ncol = *ncolx, n=*nrowx;
	double loglambda,lambda,neglik,logtau, tau;
	tau =p[ncol];
	neglik = 0;
	for(i=0; i<n; i++)
	{
	loglambda =0;
	 for(j=0; j<ncol;j++) loglambda += x[i + j*n]*p[j];
	 lambda = exp(loglambda)*ee[i];
	 //loglambda = log(lambda);
		neglik += w[i]*(lgamma(y[i]+1) + lgamma(tau) - lgamma(tau + y[i]) -
		tau * log(tau/(tau + lambda)) - y[i]* log(lambda/(tau + lambda)));
	}
	*nll = neglik;
}

void NegbinGrad(int n, double *para, double *gr, void *ex, SEXP x1){

	int i, j,ncol, N ;
	double loglambda,lambda,neglik,logtau, tau, tau_yt;
	double *exPara,*y, *x, *p, *w;
  
  exPara = (double *) ex;
  N 		= (int) exPara[0];
  ncol 	= (int)  exPara[1];
  y     = exPara + 2;
  x 		= y + N;
	w			= x + N*ncol;
	p = para;

	for(j=0; j<n;j++)  gr[j] = 0; // initialize grad
	tau =p[ncol];
	
	for(i=0; i<N; i++)
	{
	 loglambda =0;
	 for(j=0; j<ncol;j++) loglambda += x[i + j*N]*p[j];
	 lambda = exp(loglambda);//*ee[i];
	 //loglambda = log(lambda);
	 for(j=0; j<ncol;j++)
	    	 //gr[j]         += x[i + j*n]*(1-(tau + y[i])/(lambda + tau))*tau;
          gr[j]         += w[i]*(x[i + j*N]*tau*(lambda-y[i])/(lambda + tau));
          tau_yt = tau + y[i];
          gr[ncol] +=   w[i]*(-1  + (tau_yt)/(lambda + tau) + log((lambda + tau)/tau)
																+ digamma(tau) - digamma(tau_yt));
                               //+ digamma(&tau) - digamma(&tau_yt));
   }
  gr[ncol]  = tau*gr[ncol];
	//print_v(gr, 0, n);
}



double NegbinNLL(int n, double* para, void* ex, SEXP x1){

	int i, j,ncol, N ;
	double loglambda,lambda,neglik,logtau, tau, ratio;
	double *exPara,*y, *x, *p, *w;
  
  exPara = (double *) ex;
  N 		= (int) exPara[0];
  ncol 	= (int)  exPara[1];
  y     = exPara + 2;
  x 		= y + N;
	w			= x + N*ncol;
	p = para;

	tau =p[ncol];
	neglik = 0;
	for(i=0; i<N; i++)
	{
	loglambda =0;
	 for(j=0; j<ncol;j++) loglambda += x[i + j*N]*p[j];
	 lambda = exp(loglambda);//*ee[i];
	 //loglambda = log(lambda);
			neglik += w[i]*(lgamma(y[i]+1) + lgamma(tau) - lgamma(tau + y[i]) -
		tau * log(tau/(tau + lambda)) - y[i]* log(lambda/(tau + lambda)));

	}
	print_v(p, 0, n);
	Rprintf("ll: %e %e %e %e\n", neglik, ratio, lgamma(tau), lgamma(tau+y[i]));
	return(neglik);
}


void nb(double *xx, double *y, double *w, double *ee,double *x, int *Rnrowx, int *Rncolx, double *ll, int *trace){
	
	int i,l;
	int N = *Rnrowx;
	int P = *Rncolx;

	double *exPara, *y1, *x1, *w1;
	exPara = (double *) Calloc(2+N+P*N+N, double); 

	//be careful of conversion from ints to doubles
	exPara[0] = N;
	exPara[1] = P;
	y1     = exPara + 2;
	x1     = y1 + N;
	w1		 = x1 + N*P;

  int npara   = P+1; 
  int lmm     = 5; 
  int fail    = 0;
  int failA   = 0;
  int failB   = 0;
  int fncount = 0;//?
  int grcount = 0;//?
  int maxit   = 100;
  int nREPORT = 5;
			
	int nbd[npara];

	//technical parameters below:
		  double *wa, *g1;
		  int *iwa;
		  SEXP xx1;
		  PROTECT(xx1 = allocVector(REALSXP, npara));
		  //consider replacing with simple Calloc
		  wa  = (double *) S_alloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm,sizeof(double));
		  iwa = (int*) R_alloc(3*npara,sizeof(int));
		  g1 = (double *)R_alloc(npara, sizeof(double));
		
		  double gr0[npara];
		  double initPara[npara];
		  double lower[npara];
		  double upper[npara];
		  double Fmin, factr, pgtol, th0, th1, th01;
		  
		  factr = 1e7;
			pgtol = 0.0;
		  
			char msg[1023];
		  
		  for(l=0; l < npara; l++){
					nbd[l] = 2;			  	
//					initPara[l] = 0;
					lower[l] = -1000;
					upper[l] = 1000; 
			}
			//nbd[npara-1] = 1;
			lower[npara-1] = 1e-10;

	for(i=0;i<N;i++){
			for(l=0;l<P;l++){
				x1[l*N+i] = x[l*N+i];
			}
			y1[i] = y[i];
			w1[i] = w[i];
	}


      lbfgsb1(npara, lmm, xx, lower, upper, nbd, &Fmin, 
           NegbinNLL, NegbinGrad, &failA, (void*)exPara, factr, pgtol,  
           &fncount, &grcount, maxit, msg, 6, nREPORT, wa, iwa, g1,xx1);
	    if(*trace > 1) Rprintf("%d %s\n", failA, msg);
	    if (failA) {
	      if (*trace)
	        Rprintf("  i=%d, fail to fit t0 model in iteration\n", i);
	    
	      //continue;
	    }

		*ll = -Fmin;
  Free(exPara);
  UNPROTECT(1);
}

