#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "utility.h"

void calculateold(double *y, double *mu, double *pi, int *RK, int *Rn, double *a, double *forward, double *backward, double *forwardbackward, double *chsi){ 

int K=*RK; 
int n=*Rn;
int i=0; int j=0; int l=0; double sum=0;
double forwardtemp[K];
double backwardtemp[K];
double forwardbackwardtemp[K];
//initialize matrices to zero
/* not needed since passing to C from R?
for(i=0; i<n; i++;){
	for(j=0; j<K; j++;){
		forward[j*K+i]=0;
		backward[j*K+i]=0;
		forwardbackward[j*K+i]=0;
		chsi[j*K+i]=0;
	}
}
*/

//Initial observation at i=0
//forward
//for(j=0; j<K; j++){
//	forward[j*n+i] = pi[j]*dpois(y[i] , mu[j*n+i], 0);
//	backward[n*(K-j)-1-i] = 1;
//}
//Here we are logging all values now
for(j=0; j<K; j++){
	forward[j*n+i] = log(pi[j])+dpois(y[i] , mu[j*n+i], 1);
	backward[n*(K-j)-1-i] = log(1);
}





//forward and backward probabilities, 
//forward: j is current state at i, l is previous state at i-1
//backward: j is current state at i, l is future state at i-1
//for(i=1; i<n; i++){
//	for(j=0; j<K; j++){
//		for(l=0; l<K; l++){
//			forward[j*n+i] = forward[j*n+i] + forward[l*n+i-1]*a[K*l+j];
//			backward[(n-j)*K-1-i] = backward[(n-j)*K-i-1] 
//														+ backward[(n-l)*K-i]*a[K*j+l]*dpois(y[n-i] , mu[(n-l)*K-i], 0);
//		}
//		forward[j*n+i] = forward[j*n+i]*dpois(y[i] , mu[j*n+i] , 0);
//	}
//}
//Here we are logging the values to prevent underflow
for(i=1; i<n; i++){
	for(j=0; j<K; j++){
		for(l=0; l<K; l++){
			forwardtemp[l] = forward[l*n+i-1]+ log(a[K*l+j]);
			backwardtemp[K-l-1] = backward[n*(K-l)-i]+log(a[K*j+(K-l-1)])+dpois(y[n-i] , mu[n*(K-l)-i], 1);
			//if(i==1) Rprintf("%f %f %f %f %d %f %d\n", backward[n*(K-l)-i], a[K*j+(K-l-1)], dpois(y[n-i] , mu[n*(K-l)-i], 1), y[n-i], n-i, mu[n*(K-l)-i], n*(K-l)-i);
			//if(i==1) Rprintf("%d %d %f %f %f %f %d %f %d\n", j, l,forward[l*n+i-1] , a[K*l+j], +dpois(y[n-i] , mu[n*(K-l)-i], 1), y[i], i, mu[n*(K-l)-i] ,n*(K-l)-i);
		}
		//if(i==1) print_v(backwardtemp, 0, K-1);
		//if(i==1) print_v(forwardtemp, 0, K-1);
		forward[j*n+i] = logsumexp(forwardtemp, 0, K-1)+dpois(y[i] , mu[j*n+i] , 1);
		backward[n*(j+1)-1-i] = logsumexp(backwardtemp, 0, K-1);
		//Rprintf("%f\n", backward[n*(K-j)-1-i]);
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
			chsi[K*j*n+l*n+i+1] = forwardbackward[j*n+i]+log(a[j*K+l])+dpois(y[i+1] , mu[l*n+i+1], 1)+ backward[l*n+i+1]-backward[j*n+i];
			//if(i==0) Rprintf("%f %f %f %f %d %f %d\n", 0.0, a[K*j+(K-l-1)], dpois(y[i+1] , mu[l*n+i+1], 1), y[n-i], n-i,  mu[l*n+i+1], l*n+i+1);
		}
	}
}

}



void vitterbiold(double *y, double *mu, double *pi, int *RK, int *Rn, double *a, double *forward, double *backward, double *forwardbackward, double *chsi){ 

int K=*RK; 
int n=*Rn;
int i=0; int j=0; int l=0; double sum=0; double val; int idx;
double forwardtemp[K]; double forwardbackwardtemp[K];

for(j=0; j<K; j++){
	forward[j*n+i] = log(pi[j])+dpois(y[i] , mu[j*n+i], 1);
}

for(i=1; i<n; i++){
	for(j=0; j<K; j++){
		for(l=0; l<K; l++){
			forwardtemp[l] = forward[l*n+i-1]+ log(a[K*l+j]);
		}
		max(forwardtemp, 0, K-1, &val, &idx);
		forward[j*n+i] = val +dpois(y[i] , mu[j*n+i] , 1);
	}
}

for(i=n-2; i>=0; i--){
	for(j=0; j<K; j++){
			forwardbackwardtemp[j] = forward[j*n+i];
	}
	max(forwardbackwardtemp, 0, K-1, &val, &idx);
	for(j=0; j<K; j++){
			forwardbackward[idx*n+i] = 1;
			forwardbackward[(1-idx)*n+i] = 0;
	}
}


}
