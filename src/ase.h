/*
 *  ase.h
 *
 *  Created by Wei Sun on 5/25/2010.
 *
 */

#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "utility.h"


double negLogH0 (int n, double* para, void* ex);
double negLogH1 (int n, double* para, void* ex);

void negGradLogH0 (int n, double* para, double* gr, void* ex);
void negGradLogH1 (int n, double* para, double* gr, void* ex);

void ase (int* dims, double* Y1, double* Y2, double* Z, char** output, 
          double* RP_cut, int* cis_only, int* cis_distance, 
          int* eChr, int* ePos, int* mChr, int* mPos,  
          int* trace, int* succeed);
