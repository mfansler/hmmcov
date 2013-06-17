/********************************************************************
 * 
 * DNA_seq.c
 *
 * copyright (c) 2009, Wei Sun
 *
 * First written Sep 20, 2007
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * C functions for the R/lm.Tilling
 *
 * extract probe sequence, etc.
 *
 *  
 *******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <ctype.h>
#include "DNA_seq.h"
#include "utility.h"

void get_seq(char** chr_seq, int* dims, char** seq){

  int i;
  int start = dims[0]-1;
  int n = dims[1];
  
  for(i=0; i<n; i++){
    seq[0][i] = chr_seq[0][start+i];
  }
}

void get_m_seq(char** chr_seq, int* start, int* dims, char** seq){

    int i, j, m, n;
    
    m = dims[0]; // how many probes
    n = dims[1]; // length of each probe
    
    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
            seq[i][j] = chr_seq[0][start[i]-1+j];
        }
    }
}

void group_seq(char** chr_seq, int* dims, char** seq){
    
    int i, j;
    int n = dims[0]; // how many probes
    int p = dims[1]; // length of each probe
    
    for(i=0; i<n; i++){
        for(j=0; j<p; j++){
            seq[i][j] = chr_seq[0][p*i+j];
        }
    }
}


int is_lower(char** seq, int* codeR, int* dims){

    int i, j, m, n, **code;

    m = dims[0]; // how many probes
    n = dims[1]; // length of each probe
    
    reorg_int(codeR, &code, m, n);
    
    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
        	if(islower(seq[i][j])){
        		 code[i][j] = 1;
        	}
        }
    }
    
    return(1);
}

int code_AT_GC(char** seq, int* codeR, int* dims){

    int i, j, m, n, **code;

    m = dims[0]; // how many probes
    n = dims[1]; // length of each probe
    
    reorg_int(codeR, &code, m, n);
    
    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
            if(seq[i][j] == 'C' || seq[i][j] == 'c'){
                code[i][j] = 1;
            }else if(seq[i][j] == 'G' || seq[i][j] == 'g'){
                code[i][j] = 1;
            }else if(seq[i][j] == 'N' || seq[i][j] == 'n'){
                code[i][j] = -1;
            }
        }
    }
    
    return(1);
}

int code_AT_GC_var_n(char** seq, int* codeR, int* dims, int* n){

    int i, j, m, nmax, **code;

    m = dims[0]; // how many probes
    nmax = dims[1]; // length of each probe
    
    reorg_int(codeR, &code, m, nmax);
    
    for(i=0; i<m; i++){
        for(j=0; j<n[i]; j++){
            if(seq[i][j] == 'C' || seq[i][j] == 'c'){
                code[i][j] = 1;
            }else if(seq[i][j] == 'G' || seq[i][j] == 'g'){
                code[i][j] = 1;
            }else if(seq[i][j] == 'N' || seq[i][j] == 'n'){
                code[i][j] = -1;
            }
        }
    }
    
    return(1);
}


int code_ATGC(char** seq, int* codeR, int* dims){

    int i, j, m, n, **code;

    m = dims[0]; // how many probes
    n = dims[1]; // length of each probe
    
    reorg_int(codeR, &code, m, n);
    
    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
            if(seq[i][j] == 'T' || seq[i][j] == 't'){
                code[i][j] = 1;
            }else if(seq[i][j] == 'C' || seq[i][j] == 'c'){
                code[i][j] = 2;
            }else if(seq[i][j] == 'G' || seq[i][j] == 'g'){
                code[i][j] = 3;
            }
        }
    }
    
    return(1);
}


