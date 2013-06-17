/**********************************************************************
 *
 * baseQ.c
 *
 * summarize the quality of each base
 *
 * copyright (c) 2009, Wei Sun, UNC-CH
 *
 * last modified May 15, 2009
 * first written May 15, 2009
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 **********************************************************************/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>

#define MAX_LEN 1023

/**********************************************************************
 *
 * generate counts from the ouptut of maq mapview
 *
 **********************************************************************/

/**
 *
 one line maq mapview output looks like:
 [1-4]   HWI-EAS156:1:82:1166:1191#0/1   chr10   3000442 +       
 [5-9]   0       0       0       0       0       
 [10-14] 0       0       2       27      35       
 [15]    AAGtTCCCAGAAAAGCTGTTTTCCTGCTACTCTCT     
 [16]    @AB;AAAABB@>A@BABB@ABB@@?BA@@=B>A>@

 Each line of maq mapview output consists of
 [1] read name
 [2] chromosome
 [3] position
 [4] strand
 [5] insert size from the outer coorniates of a pair
 [6] paired flag
 [7] mapping quality
 [8] single-end mapping quality
 [9] alternative mapping quality
 [10] number of mismatches of the best hit
 [11] sum of qualities of mismatched bases of the best hit
 [12] number of 0-mismatch hits of the first 24bp
 [13] number of 1-mismatch hits of the first 24bp on the reference
 [14] length of the read
 [15] read sequence
 [16] its quality
 */

void baseQ(char **Rinput, char **Routput, int* RL, int* skip,
 double *Rm1, double *Rm2, double* Rmin, double *Rmax)
{
  int L = *RL, out, i, j, k;
  double tmp;
  char str_buf[MAX_LEN];
  char *input, *output;
	char *ptr= NULL;
	char *delim = "\t";
	char qua[255];
  FILE *FI, *FO;

  input   = Rinput[0];
  output  = Routput[0];

  if(strcmp(output, "")==0){
    out = 0;
  }else{
    out = 1;
  }

  FI = fopen(input, "r");
  assert(FI);
  
  if(out){
    FO = fopen(output, "w");
    assert(FO);
  }
  
  for(i=0; i < *skip; i++){
    if(fgets(str_buf, MAX_LEN, FI) == NULL){
      error("thre are only %d lines in file %s\n", i, input);
    }
  }

  /**
   * read in the first line
   * obtain the initial value for min and max
   */
  if(fgets(str_buf, MAX_LEN, FI) == NULL){
    error("file %s is empty after skipping %d lines\n", input, skip);
  }
  i = 0;
  if ((ptr = strtok(str_buf, delim)) != NULL) {
    do {
      i++;
      if(i==16){
        strcpy(qua, ptr);
      }
    } while ((ptr = strtok(NULL, delim)) != NULL);
  }else{
    error("%s is not tab-delimated\n", input);
  }
  
  for(j=0; j<L; j++){      
    tmp = (double)qua[j];
    if(out){
      fprintf(FO, "%d ", (int)tmp);
    }
    Rm1[j] = tmp;
    Rm2[j] = tmp*tmp;
    Rmax[j] = tmp;
    Rmin[j] = tmp;
  }
  if(out){
    fprintf(FO, "\n");
  }

  /**
   * read in the the other lines
   */

  k = 1;
  while(fgets(str_buf, MAX_LEN, FI) != NULL){
    k++;
    i = 0;
    if ((ptr = strtok(str_buf, delim)) != NULL) {
      do {
        i++;
        if(i==16){
          strcpy(qua, ptr);
        }
      } while ((ptr = strtok(NULL, delim)) != NULL);
    }else{
      error("%s is not tab-delimated\n", input);
    }
    
    for(j=0; j<L; j++){      
      tmp = (double)qua[j];
      if(out){
        fprintf(FO, "%d ", (int)tmp);
      }
      Rm1[j] += tmp;
      Rm2[j] += tmp*tmp;
      
      if(Rmax[j] < tmp){ Rmax[j] = tmp; }
      if(Rmin[j] > tmp){ Rmin[j] = tmp; } 
    }
    if(out){
      fprintf(FO, "\n");
    }
  }
  
  for(j=0; j<L; j++){
    Rm1[j] /= (double)k;
    Rm2[j] /= (double)k;
  }
  
  fclose(FI);
  if(out) fclose(FO);
  
}
