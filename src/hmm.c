#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "utility.h"

void calculate(double *y, double *mu, double *pi, int *RK, int *Rn, double *a, double *forward, double *backward, double *forwardbackward, double *chsi){ 

int K=*RK; 
int n=*Rn;
int i=0; int j=0; int l=0; int idx; double sum=0;
double forwardtemp[K];
double backwardtemp[K];
double forwardbackwardtemp[K];
//initialize matrices to zero
/* not needed since passing to C from R?*/
for(i=0; i<n; i++){
	for(j=0; j<K; j++){
		forward[j*K+i]=0;
		backward[j*K+i]=0;
		forwardbackward[j*K+i]=0;
		for(l=0; l<K; l++) chsi[K*j*n+l*n+i]=0;
	}
}

	for(j=0; j<K; j++){
		forwardtemp[j]=0;
		backwardtemp[j]=0;
		forwardbackwardtemp[j]=0;
	}

//Initial observation at i=0
//forward
//for(j=0; j<K; j++){
//	forward[j*n+i] = pi[j]*dpois(y[i] , mu[j*n+i], 0);
//	backward[n*(K-j)-1-i] = 1;
//}
//Here we are logging all values now
i=0;
for(j=0; j<K; j++){
	forward[j*n+i] = log(pi[j])/*pi[j]*/+dpois(y[i] , mu[j*n+j*n+i], 1);
	backward[n*(K-j)-1-i] = log(1);
}


for(i=1; i<n; i++){
	for(j=0; j<K; j++){
		for(l=0; l<K; l++){
			forwardtemp[l] = forward[l*n+i-1] + log(a[K*l+j])/*a[K*l+j]*/ + dpois(y[i] , mu[n*K*j+l*n+i] , 1);
			backwardtemp[l] = backward[n*(K-1) + l*n-i]+log(a[K*j+l]) /*a[K*j+l]*/ +dpois(y[n-i] , mu[n*(K-1)+n*(K-1)*j+n*K*l-i], 1);
		}
		forward[j*n+i] = logsumexp(forwardtemp, 0, K-1);
		backward[n*(K-1) + j*n-i-1] = logsumexp(backwardtemp, 0, K-1);
	}
}
 
//forward-backward probabilitues
for(i=0; i<n; i++){
	for(j=0; j<K; j++){
			forwardbackwardtemp[j] = forward[j*n+i]+backward[j*n+i];
	}
	sum = logsumexp(forwardbackwardtemp, 0, K-1);
	//print_v(forwardbackwardtemp, 0, K-1);
	//Rprintf("%f\n",sum);
	for(j=0; j<K; j++){
			forwardbackward[j*n+i] = forwardbackwardtemp[j]-sum;
	}
}

//chsi (p(Qt=i, Qt+q=j, O | lambda))
for(i=0; i<n-1; i++){
	for(j=0; j<K; j++){
		for(l=0; l<K; l++){
			chsi[K*j*n+l*n+i] = forwardbackward[j*n+i]+log(a[j*K+l]) /*a[j*K+l]*/ + dpois(y[i+1] , mu[n*(K-1)*j+n*K*l+i+1], 1)+ backward[l*n+i+1]-backward[j*n+i];
		}
	}
}
/*
for(i=0; i<n-1; i++){
	for(j=0; j<K; j++){
		for(l=0; l<K; l++){
			chsi[K*j*n+l*n+i] = forwardbackward[j*n+i]+/*log(a[j*K+l]) a[j*K+l] + dpois(y[i+1] , mu[n*(K-1)*j+n*K*l+i+1], 1)+ backward[l*n+i+1]-backward[j*n+i];
		}
	}
}
*/

}


void vitterbi(double *y, double *mu, double *pi, int *RK, int *Rn, double *a, double *forward, double *backward, double *forwardbackward, double *chsi){ 

int K=*RK; 
int n=*Rn;
int i=0; int j=0; int l=0; double sum=0; double val=0; int idx=0;
double forwardtemp[K];
double backwardtemp[K];
double forwardbackwardtemp[K];

for(j=0; j<K; j++){
	forward[j*n+i] = /*log(pi[j])*/ pi[j] +dpois(y[i] , mu[j*n+j*n+i], 1);
}


for(i=1; i<n; i++){
	for(j=0; j<K; j++){
		for(l=0; l<K; l++){
			forwardtemp[l] = forward[l*n+i-1] + /*log(a[K*l+j])*/ a[K*l+j] + dpois(y[i] , mu[n*K*j+l*n+i] , 1);
		}
		max(forwardtemp, 0, K-1, &val, &idx);
		forward[j*n+i] = val ;
	}
}


for(i=n-2; i>=0; i--){
	for(j=0; j<K; j++){
			forwardbackwardtemp[j] = forward[j*n+i];
	}
	max(forwardbackwardtemp, 0, K-1, &val, &idx);
	for(j=0; j<K; j++){
			forwardbackward[idx*n+i] = 1;
	}
}

}

