/**********************************************************************
 * 
 * ss_hmm.c
 *
 * copyright (c) 2006, Wei Sun, UNC
 *
 * last modified Dec 10, 2009
 * first written Summer, 2006
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/ss.hmm package
 *
 * Segmental Semi-Markov Hidden Markov Model
 *
 *  
 **********************************************************************/


/**********************************************************************
 *   HMM: 4 states: 
 *
 * state 0
 *   flat line long (nucleosome)
 *       emiss: dnorm(0, sigam1)
 *
 * state 1
 *   flat line short (nucleosome free)
 *       emiss: dnorm(0, sigam2)
 *
 * state 2
 *   line of negative slope (followed by flat line short)
 *       emiss: dnorm(0, sigam3)
 *
 * state 3
 *   line of positive slope (followed by flat line long)
 *       emiss: dnorm(0, sigam4)
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include "utility.h"
#include "ss_hmm.h"

/**********************************************************************
 * 
 * get_end
 *
 * get the fitted value for the end of one segment, given a fix point. 
 *
 * The fix point could be the end of previous segment, or the
 * the mean value of x and p
 *
 * this function works only for the four states: 
 * 0: long flat line
 * 1: short flat line
 * 2: negative slope line
 * 3: positive slope line
 *
 * NOTICE: WE ONLY CALCULATE THE ENDS IF IT IS POSSIBLE 
 * IMPOSSIBLE fitted value at the ENDES ARE SET AS 0.0 AS DEFAULT
 * THOSE ENDS ARE NOT SUPPOSED TO BE USED IN OTHER FUNCTIONS
 * 
 * ends_idx indicate those ends of impossible segments:
 * segment (t1, t2, i) is impossible if emission probability = -inf
 * duration probability is not considered here, it is accounted in 
 * function get_ends
 *
 * ends: a vector of four elements, which stores the fitted values of 
 *       of the four states
 *
 **********************************************************************/

void get_end(double* x, double* p, double p_end, int nl, int nh, 
             double x_fix, double p_fix, int state, double* ends, 
             int* ends_idx, double bLimit)
{
  int k, j;
  int n = nh - nl + 1;
  double *nx, *np, b;
    
	nx = (double *) calloc(n, sizeof(double));
	np = (double *) calloc(n, sizeof(double));
    
  if(state<0){ // calculate ends of all states

    ends[0] = ends[1] = x_fix;
    ends_idx[0] = ends_idx[1] = 1;
    
    j = nl;
    for(k=0; k<n; k++){
      nx[k] = x[j] - x_fix;
      np[k] = p[j] - p_fix; 
      j++;
    }
    b = slope(nx, np, n);
    
    if(fabs(b) <= bLimit){
      ends[2] = 0.0;
      ends[3] = 0.0;
      ends_idx[2] = 0;
      ends_idx[3] = 0;
    }else if(b < 0){
      ends[2]     = x_fix + b*(p_end - p_fix);
      ends[3]     = 0.0;
      ends_idx[2] = 1;
      ends_idx[3] = 0;
    }else{
      ends[2]     = 0.0;
      ends[3]     = x_fix + b*(p_end - p_fix);
      ends_idx[2] = 0;
      ends_idx[3] = 1;
    }
    
  }else{ // calculate end of one state

    if(state==0){
      ends[0]     = x_fix;
      ends_idx[0] = 1;
    }else if(state==1){
      ends[1]     = x_fix;
      ends_idx[1] = 1;
    }else {
      j = nl;
      for(k=0; k<n; k++){
        nx[k] = x[j] - x_fix;
        np[k] = p[j] - p_fix;
        j++;
      }
      b = slope(nx, np, n);
        
      if(state==2){
        if(b < -bLimit){
          ends[2]     = x_fix + b*(p_end - p_fix);
          ends_idx[2] = 1;
        }else{
          ends[2]     = 0.0;
          ends_idx[2] = 0;
        }
      }else{
        if(b > bLimit){
          ends[3]     = x_fix + b*(p_end - p_fix);
          ends_idx[3] = 1;
        }else{
          ends[3]     = 0.0;
          ends_idx[3] = 0;
        }
      }
    }
  }
  
  free(nx);
  free(np);
}

/**********************************************************************
 * 
 * get_ends
 *
 * get the ends of all segments
 * 
 * this function works only for the four states: 
 * 0: long flat line
 * 1: short flat line
 * 2: negative slope line
 * 3: positive slope line
 *
 * NOTICE: WE ONLY CALCULATE THE ENDS IF IT IS POSSIBLE 
 * IMPOSSIBLE ENDES ARE SET AS 0.0 AS DEFAULT
 * THOSE ENDS ARE NOT SUPPOSED TO BE USED IN OTHER FUNCTIONS
 *
 * xR_end_idx indicate those ends of impossible segments
 * segment (t1, t2, i) is impossible if 
 * (1) emission probability = -inf, or 
 * (2) duration probability = -inf
 **********************************************************************/

void get_ends(double* x, double* p, int* stp_end, int* dims, double* duR, 
              double* xR_end, int* xR_end_idx, double bLimit)
{
  int i, j, j_end, d, s, n, m, D, Z, L, nl, nh, *endIdx, **x_end_idx;
  double p_end, xf, pf, *ends, **x_end, **du;
  
  n = dims[0]; // number of observations
  m = dims[1]; // number of states
  D = dims[2]; // longest duration
  Z = dims[3]; // number of bins
  L = dims[4]; // length of each step
  
  reorg(duR, &du, D, m);
  reorg(xR_end, &x_end, m*Z, Z);
  reorg_int(xR_end_idx, &x_end_idx, m*Z, Z);

  /*
  Rprintf("duration probability matrix:\n");
  Rprint_me(du, 0, 10, 0, m-1);
  
  Rprintf("x_end:\n");
  Rprint_me(x_end, 0, 10, 0, 3);
  
  Rprintf("x_end_idx:\n");
  Rprint_mi(x_end_idx, 0, 10, 0, 3);

  Rprintf("stp_end:\n");
  Rprint_vi(stp_end, 0, 20);
  */
  
  /* ends at one time point */
	ends   = (double *)calloc(m, sizeof(double));
	endIdx = (int *)calloc(m, sizeof(int));
    
  i = 0;
  j_end = min_int(D-1, Z-1);
  
  for(j=0; j<=j_end; j++){
    nl = 0;
    nh = stp_end[j];
    d  = j - i; // d = actual duration - 1
    
    xf = mean(x, nl, nh);
    pf = mean(p, nl, nh);
    p_end = (j+1)*L;

    get_end(x, p, p_end, nl, nh, xf, pf, -1, ends, endIdx, bLimit);
    
    for(s=0; s<m; s++){
      if(du[d][s]>0){ // if duration prob > 0
        x_end[s*Z][j] = ends[s];
        x_end_idx[s*Z][j] = endIdx[s];
      }
    }
  }
  
  for(i=1; i<Z; i++){
    j_end = min_int(i+D-1, Z-1);
    for(j=i; j<=j_end; j++){
      nl = stp_end[i-1]+1;
      nh = stp_end[j];
      d  = j - i; // d = actual duration - 1
      xf = mean(x, nl, nh);
      pf = mean(p, nl, nh);
      p_end = (j+1)*L;
      get_end(x, p, p_end, nl, nh, xf, pf, -1, ends, endIdx, bLimit);
      for(s=0; s<m; s++){
        if(du[d][s]>0){ // if duration prob > 0
          x_end[s*Z+i][j] = ends[s];
          x_end_idx[s*Z+i][j] = endIdx[s];
        }
      }
    }
  }
  
  free(ends);
  free(endIdx);
}

/**********************************************************************
 * 
 * emiss
 *
 * (log) emission probability, with fix point;  
 *
 * the fix point could be the end of the previous segment or 
 * the mean value of x and p
 *
 **********************************************************************/

int emiss(double* x, double* p, int nl, int nh, double x_fix, double p_fix, 
          int* state, double* ev, int penalty, double bLimit)
{
  double sd_x, sd_r, tmp, b, s, pena;
  int i, j;
  int n = nh - nl + 1;
  double *nx, *np, *resid;
      
  if(n<3){
    error("\n number of observations < 3 in function emiss \n");
  }
  
	nx = (double *) calloc(n, sizeof(double));
	np = (double *) calloc(n, sizeof(double));
	resid = (double *) calloc(n, sizeof(double));
    
  if(state[0] || state[1]){
    /* ----------------------------------------------------------
    * state 0 or 1, flat lines
    * ----------------------------------------------------------*/
    sd_x = sd(x, nl, nh, x_fix);
    
    tmp  = 0.0;
    for(i=nl; i<=nh; i++){
      tmp += dnorm(x[i], x_fix, sd_x, 1); //dnorm(, log=TRUE)
    }

    if(state[0]){ ev[0] = tmp; }
    if(state[1]){ ev[1] = tmp; }
  }
  
  if(state[2] || state[3]){
    /* ----------------------------------------------------------
    * state 2 or 3, line with slope
    *
    * x = a + b*p             (1)
    * x.bar = a + b*p.bar     (2)
    * x.fix = a + b*p.fix     (3)
    * 
    * (1) - (2) => x-x.bar = b*(p - p.bar)
    * (1) - (3) => x-x.fix = b*(p - p.fix)
    * ----------------------------------------------------------*/
    
    j = 0;
    for(i=nl; i<=nh; i++){
      nx[j] = x[i] - x_fix;
      np[j] = p[i] - p_fix;
      j++;
    }
    
    b = slope(nx, np, n);
    
    if(fabs(b) <=  bLimit){ b = 0.0; }
    
    if(penalty==0){ /* no penalty */
      pena = 0.0;
    }else if(penalty==2){ /* BIC */
      pena = 0.5*log(n);
    }else if(penalty==1){ /* AIC */
      pena = 1;
    }else{
      error("invalid penalty\n");
    }
    
    /*
     * 0: long flat line
     * 1: short flat line
     * 2: negative slope line
     * 3: positive slope line
     */
    
    if(!state[3]){
      // if only calculate emission probability of state 2
      if(b >= -bLimit){
        ev[2] = log(0.0);
      }else{
        for(i=0; i<n; i++){ resid[i] = nx[i] - b*np[i]; }
        sd_r  = sd(resid, 0, n-1, 0.0);
        ev[2] = logL(resid, n, 0.0, sd_r) - pena;
      }
    }
    
    if(!state[2]){
      // if only calculate emission probability of state 3
      if(b <= bLimit){
        ev[3] = log(0.0);
      }else{
        for(i=0; i<n; i++){ resid[i] = nx[i] - b*np[i]; }
        sd_r  = sd(resid, 0, n-1, 0.0);
        ev[3] = logL(resid, n, 0.0, sd_r) - pena;
      }
    }
    
    if(state[2] && state[3]){
      // calculate emission probability of both state 2 and 3
      for(i=0; i<n; i++){ resid[i] = nx[i] - b*np[i]; }
      sd_r = sd(resid, 0, n-1, 0.0);
      
      if(fabs(b) <= bLimit){
        ev[2] = log(0.0);
        ev[3] = log(0.0);
      }else if(b > 0){
        ev[2] = log(0.0);
        ev[3] = logL(resid, n, 0.0, sd_r) - pena;
      }else{
        ev[3] = log(0.0);
        ev[2] = logL(resid, n, 0.0, sd_r) - pena;
      }
    }
  }
  
  free(nx);
  free(np);
  free(resid);
  
  return(1);
}

/**********************************************************************
 * 
 * viterbi
 *
 * find the best path
 *
 **********************************************************************/


int viterbi(double* x, double* p, int* stp_end, int* dims, double* path, 
            double* segmentR, int* n_seg, double* log_lik, double* aR, 
            double* duR, int* diR, double* pi, int* penalty, double* bLimit)
{
  int i, j, t, tp, k, d, d_end, t1, s1, t2, s2, d_idx, s_idx;
  int nl, nh, which_i, which_d, seg, p_end;
  int *state, *endIdx, **prev, **dura, **di;
  double xf, pf, max_lik, lik_t, d_lik, *lik, xbegin;
  double **xend, **logp, *ends, *ev, **em;
  double **segment, **a, **du, *log_pi, **log_a, **log_du;
  
  int n, m, D, Z, L;
  n = dims[0]; // number of observations
  m = dims[1]; // number of states
  D = dims[2]; // longest duration
  Z = dims[3]; // number of steps
  L = dims[4]; // length of each step
  
  reorg(segmentR, &segment, Z+1, 3);
  reorg(aR, &a, m, m);
  reorg(duR, &du, D, m);
  reorg_int(diR, &di, D, m);

  /*
  Rprintf("Z=%d, penalty=%d, bLimit=%f\n", Z, *penalty, *bLimit);
  Rprintf("transition probability matrix\n");
  Rprint_me(a, 0, m-1, 0, m-1);
  Rprintf("duration probability \n");
  Rprint_me(du, 0, 20, 0, m-1);

  Rprintf("duration probability index\n");
  Rprint_mi(di, 0, 20, 0, m-1);
  */

  /* state */
	state  = (int *)calloc(m, sizeof(int));
	
  /* ends at one time point */
	ends   = (double *)calloc(m, sizeof(double));
	endIdx = (int *)calloc(m, sizeof(int));

  /* likelihood for different duration */
	lik  = (double *)calloc(Z, sizeof(double));

  /* emission value */
	ev = (double *)calloc(m, sizeof(double));
  
	em = (double**) malloc(m*Z*sizeof(double*));
  em[0] = (double*) calloc(m*Z*m, sizeof(double));
  for(i=1; i<m*Z; i++){
    em[i] = em[0] + i*m;
  }

  /* dura[t][i]: duration of state i which ends at step t */
  dura = (int **)malloc(Z*sizeof(int*));
	dura[0] = (int*) calloc(Z*m, sizeof(int));
	for(i=1; i<Z; i++){
    dura[i] = dura[0] + i*m;
	}
	
  /* prev[t][i]: previous state of state i which ends at step t */
  prev = (int **)malloc(Z*sizeof(int*));
	prev[0] = (int*) calloc(Z*m, sizeof(int));
	for(i=1; i<Z; i++){
    prev[i] = prev[0] + i*m;
	}
	
  /* xend[t][i]: the expected ending signal of state i which ends at step t */
  xend = (double **)malloc(Z*sizeof(double*));
	xend[0] = (double*) calloc(Z*m, sizeof(double));
	for(i=1; i<Z; i++){
    xend[i] = xend[0] + i*m;
	}
	
  /* logp[t][i]: log likelihood if the path ends at state i, step t */
  logp = (double **)malloc(Z*sizeof(double*));
	logp[0] = (double*) calloc(Z*m, sizeof(double));
	for(i=1; i<Z; i++){
    logp[i] = logp[0] + i*m;
	}

  /* log probability */
  log_pi = (double *)calloc(m, sizeof(double));

	log_a  = (double**) malloc(m*sizeof(double*));
	log_a[0] = (double*) calloc(m*m, sizeof(double));
	for(i=1; i<m; i++){
    log_a[i] = log_a[0] + i*m;
	}
	
	log_du = (double**) malloc(D*sizeof(double*));
	log_du[0] = (double*) calloc(D*m, sizeof(double));
	for(i=1; i<D; i++){
    log_du[i] = log_du[0] + i*m;
	}

  /* ---------------------------------------------------
   * log transformation of the probability 
   * --------------------------------------------------- */
  for(i=0; i<m; i++){
    log_pi[i] = log(pi[i]);
    for(j=0; j<m; j++){
      log_a[i][j] = log(a[i][j]);
    }
  }
  
  for(i=0; i<D; i++){
    for(j=0; j<m; j++){
      log_du[i][j] = log(du[i][j]);
    }
  }
    
  /* ---------------------------------------------------
   * Initialization 
   * --------------------------------------------------- */
  
  /**
   * the first stp_end[0] data points belong to the  
   * first bin. We make prediction of the ends at 
   * position L, the end of the first bin
   */
  
  nl = 0; 
  nh = stp_end[0];
  p_end = L;
  
  /* state[i] is inidicator of which state to to calculate emission prob */
  for(i=0; i<m; i++){ state[i] = 1;}
  
  xf = mean(x, nl, nh);
  pf = mean(p, nl, nh);
  emiss(x, p, nl, nh, xf, pf, state, ev, *penalty, *bLimit);
  get_end(x, p, p_end, nl, nh, xf, pf, -1, ends, endIdx, *bLimit);
  
  for(i=0; i<m; i++){
    logp[0][i] = log_pi[i] + log_du[0][i] + ev[i];
    /*
     * we do not need to use endIdx to indicate whether there is 
     * a valid end, dura[t][i] and ev can do this
     */
    
    if(is_infinite(logp[0][i]) == - 1){
      dura[0][i] = -1;
      xend[0][i] = 0.0;
    }else{
      dura[0][i] = 0; // dura[t][j] = actual duration - 1
      xend[0][i] = ends[i];
    }
    prev[0][i] = -1;
  }
  
  /* ---------------------------------------------------
   * Recursion 
   * --------------------------------------------------- */
  for(t=1; t<Z; t++){
    p_end = (t+1)*L;
    d_end = min_int(D-1, t-1); // d = actual duration - 1
    nh    = stp_end[t];
    
    /* 
     * calculate emission probability 
     * save emission prob for different i and d
     */

    for(i=0; i<m; i++){
      // state j after state i
      for(d=0; d<=d_end; d++){
        for(j=0; j<m; j++){
          if((a[i][j]>0) && (du[d][j]>0) ){ state[j] = 1; }
          else{ state[j] = 0; }
        }
        tp = t - (d+1); // end step of state i
        nl = stp_end[tp] + 1;
        pf = L*(tp+1);
        xf = xend[tp][i];  // fixed value of x
        /* em[i*Z+d][j] = emiss prob of emiss j, from state i to j */
        emiss(x, p, nl, nh, xf, pf, state, em[i*Z+d], *penalty, *bLimit);
      }
    }

    /* calculate likelihood of each path */
    for(j=0; j<m; j++){
        
      max_lik = log(0.0);
      which_i = -1;
      which_d = -1;
      
      for(i=0; i<m; i++){
        if(a[i][j] > 0){
          for(d=0; d<=d_end; d++){
            /* sometime the em[i*Z+d][j] is not calculated
             * because either a[i][j] = 0 or du[d][j] = 0
             */
            tp = t - (d+1); // end step of state i
            lik[d] = logp[tp][i] + log_a[i][j] + log_du[d][j] 
                        + em[i*Z+d][j];
          }
          
          max(lik, 0, d_end, &d_lik, &d_idx);
          
          if(d_lik > max_lik){ 
            max_lik = d_lik;
            which_d = d_idx; 
            which_i = i;
          }
        }
          
      }
      
      logp[t][j] = max_lik;
      if(is_infinite(max_lik) == - 1){
        dura[t][j] = -1;
        prev[t][j] = -1;
      }else{
        dura[t][j] = which_d;
        prev[t][j] = which_i;
      }
        
    }
    
    /*
     * if duration=t, there is no previous state i
     */
    if(D > t){ // t = actual step - 1, D>t <=> D>=t+1
      nl = 0;
      nh = stp_end[t];
      xf = mean(x, nl, nh);
      pf = mean(p, nl, nh);
      for(k=0; k<m; k++){state[k]=1;}
      emiss(x, p, nl, nh, xf, pf, state, ev, *penalty, *bLimit);
      for(j=0; j<m; j++){
        // in log_du[t][j], t=actual duration - 1
        lik_t = log_pi[j] + log_du[t][j] + ev[j];
        if(lik_t > logp[t][j]){
          logp[t][j] = lik_t;
          dura[t][j] = t; // dura[t][j] = actual duration - 1
          prev[t][j] = -1;
        }
      }
    }
    
    /*
     * find the fix ponit
     */
    for(j=0; j<m; j++){
      d = dura[t][j];
      if(d>=0){
        nh = stp_end[t];
        /*
         * since d>=0, it is guaranteed that there is valid ends
         * so we do not need to use endIdx here
         */
        if(d==t){ 
          nl = 0; 
          xf = mean(x, nl, nh);
          pf = mean(p, nl, nh);
          get_end(x, p, p_end, nl, nh, xf, pf, j, ends, endIdx, *bLimit);
          xend[t][j] = ends[j];
        }else {
          tp = t - (d+1); // end step of the state before state j
          nl = stp_end[tp]+1;
          pf = (tp+1)*L;
          i  = prev[t][j]; // the state before state j
          xf = xend[tp][i];  // fixed value of x
          get_end(x, p, p_end, nl, nh, xf, pf, j, ends, endIdx, *bLimit);
          xend[t][j] = ends[j];
        }
      }
    }
      
  }

  seg = Z+1;
  max(logp[Z-1], 0, m-1, log_lik, &s_idx);
  if(is_infinite(*log_lik) == - 1){
    error("no path found in viterbi\n");
  }else{
    //Rprintf("log_lik = %f\n", *log_lik);
    t2 = Z-1;
    s2 = s_idx;
    seg--;
    segment[seg][0] = (t2+1)*L;
    segment[seg][1] = xend[t2][s2];
    segment[seg][2] = s2;
    d  = dura[t2][s2];
    s1 = prev[t2][s2];
    t1 = t2 - d; // d = actual duration - 1, t1 = start step of the segment
    for(t=t1; t<=t2; t++){
      path[t] = s2;
    }
  }
  
  while(t1>0){
    t2 = t1-1; // the end of one segment
    s2 = s1;
    seg--;
    segment[seg][0] = (t2+1)*L;
    segment[seg][1] = xend[t2][s2];
    segment[seg][2] = s2;
    d  = dura[t2][s2];
    s1 = prev[t2][s2];
    if(d<0){
      warning("incomplete path found\n");
    }
    t1 = t2 - d; // d = actual duration - 1, t1 = start step of the segment
    for(t=t1; t<=t2; t++){
      path[t] = s2;
    }
  }
  
  /* add the beginning point of the fitted line */
  nl = 0;
  nh = stp_end[t2];
  pf = segment[seg][0];
  xf = segment[seg][1];
  p_end = 1;
  j  = path[0];
  get_end(x, p, p_end, nl, nh, xf, pf, j, ends, endIdx, *bLimit);
  seg--;
  segment[seg][0] = 1;
  segment[seg][1] = ends[j];
  segment[seg][2] = j;
  
  *n_seg = Z+1 - seg;
  
  free(state);
  free(ends);
  free(endIdx);
  free(lik);
  free(ev);
  free(em);
  free(prev);
  free(dura);
  free(xend);
  free(logp);
  free(log_pi);
  free(log_a);
  free(log_du);
  
  return(1);
}

