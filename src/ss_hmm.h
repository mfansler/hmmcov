/**********************************************************************
 * 
 * hmm_linear.h
 *
 * copyright (c) 2006, Wei Sun, UCLA
 *
 * last modified Sep 27, 2006
 * first written Aug 29, 2006
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/ss.hmm package
 *
 * Segmental Semi-Markov Hidden Markov Model
 *
 **********************************************************************/
 
void get_end(double* x, double* p, double p_end, int nl, int nh, double x_fix, double p_fix, int state, double* ends, int* ends_idx, double bLimit);

void get_ends(double* x, double* p, int* stp_end, int* dims, double* duR, double* xR_end, int* xR_end_idx, double bLimit);

int emiss(double* x, double* p, int nl, int nh, double x_fix, double p_fix, int* state, double* ev, int penalty, double bLimit);

int forward(double* x, double* p, int* stp_end, int* dims, double* fR, double* aR, double* duR, int* diR, double* xR_end, int* xR_end_idx, double* pi, int* penalty, double* bLimit);

int backward(double* x, double* p, int* stp_end, int* dims, double* bR, double* aR, double* duR, int* diR, double* xR_end, int* xR_end_idx, int* penalty, double* bLimit);

int post_prob(double* fR, int* fR_inf, double* bR, int* bR_inf, double* postPR, int* dims);

int viterbi(double* x, double* p, int* stp_end, int* dims, double* path, double* segmentR, int* n_seg, double* log_lik, double* aR, double* duR, int* diR, double* pi, int* penalty, double* bLimit);
