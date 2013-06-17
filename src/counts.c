/**********************************************************************
 *
 * counts.c
 *
 * allele specfic sequecing study
 *
 * copyright (c) 2009, Wei Sun, UNC-CH
 *
 * last modified May 01, 2009
 * first written May 01, 2009
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>

#define MAX_LEN 1025

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

void counts(char **Rinput, char **Routput, int *Rm, int *RQ, int *Rq, 
 int *Rnl0, int *Rnu0, int *Rn1, int *skip)
{
/**
 *
  m INT 	Maximum number of mismatches allowed for a read to be used [4]
  Q INT 	Maximum allowed number of quality values of mismatches [50]
  q INT 	Minimum mapping quality allowed for a read to be used [0] 
  nl0 INT Minimum allowed number of 0-mismatch hits of the first 24bp [1]
  nu0 INT Maximum allowed number of 0-mismatch hits of the first 24bp [1]
  n1  INT Maximum allowed number of 1-mismatch hits of the first 24bp [5]
 */
  int m = *Rm;
  int Q = *RQ;
  int q = *Rq;
  int nl0 = *Rnl0;
  int nu0 = *Rnu0;
  int n1  = *Rn1;
  int i, maq, nm, sqm, k5, k3, nm0, nm1;
  long int pos, pos1;
  FILE *FI, *FO;
  char str_buf[MAX_LEN], line[MAX_LEN];
  char chr[128], chr1[128], strand[2];
	char *delim = "\t";
	char *ptr= NULL;
  char *input, *output;
  input  = Rinput[0];
  output = Routput[0];
  
  // Rprintf("input is %s, and output is %s\n", input, output);
  
  FI = fopen(input, "r");
  if(FI == NULL){
    error("cannot open file %s\n", input);
  }
  
  FO = fopen(output, "w");
  if(FO == NULL){
    error("cannot open file %s\n", output);
  }
  
  for(i=0; i < *skip; i++){
    if(fgets(str_buf, MAX_LEN, FI) == NULL){
      error("thre are only %d lines in file %s\n", i, input);
    }
  }

  /**
   * read in the first line
   * obtain the initial chromosome (chr1) and position (pos1)
   */
  k5 = k3 = 0;

  if(fgets(str_buf, MAX_LEN, FI) == NULL){
    error("file %s is empty after skipping %d lines\n", input, skip);
  }
    
  i = 0;
  if ((ptr = strtok(str_buf, delim)) != NULL) {
    do {
      i++;
      if(i==2){
        strcpy(chr1, ptr);
      }else if(i==3){
        pos1 = atol(ptr);
      }else if(i==4){
        strcpy(strand, ptr);
      }else if(i==7){
        maq = atoi(ptr);
      }else if(i==10){
        nm  = atoi(ptr);
      }else if(i==11){
        sqm = atoi(ptr);
      }else if(i==12){
        nm0 = atoi(ptr);
      }else if(i==13){
        nm1 = atoi(ptr);
      }
    } while ((ptr = strtok(NULL, delim)) != NULL);
  }else{
    error("%s is not tab-delimated\n", input);
  }

  //Rprintf("chr: %s, pos: %d, strand: %s, chr1, pos1, strand);
  //Rprintf("score: %d, nmismatch: %d, score of mismatch: %d\n", maq, nm, sqm);
  
  if(maq >= q && nm <= m && sqm <= Q && nm0 <= nu0 && nm0 >= nl0 && nm1 <= n1){
    if(strcmp(strand, "+")==0){
      k5++;
    }else if(strcmp(strand, "-")==0){
      k3++;
    }else{
      error("unrecognized strand %s \n", strand);
    }
  }

  while(fgets(str_buf, MAX_LEN, FI) != NULL){
    i = 0;
    if ((ptr = strtok(str_buf, delim)) != NULL) {
      do {
        i++;
        if(i==2){
          strcpy(chr, ptr);
        }else if(i==3){
          pos = atol(ptr);
        }else if(i==4){
          strcpy(strand, ptr);
        }else if(i==7){
          maq = atoi(ptr);
        }else if(i==10){
          nm  = atoi(ptr);
        }else if(i==11){
          sqm = atoi(ptr);
        }else if(i==12){
          nm0 = atoi(ptr);
        }else if(i==13){
          nm1 = atoi(ptr);
        }
      } while ((ptr = strtok(NULL, delim)) != NULL);
    }else{
      error("%s is not tab-delimated\n", input);
    }
    
    if(maq >= q && nm <= m && sqm <= Q && nm0 <= nu0 && nm0 >= nl0 && nm1 <= n1)
    {
      /* check whether we are moving to the next posiction */
      if(strcmp(chr, chr1) == 0 && pos == pos1){
        if(strcmp(strand, "+")==0){
          k5++;
        }else if(strcmp(strand, "-")==0){
          k3++;
        }else{
          error("unrecognized strand %s \n", strand);
        }
      }else{
        if(strcmp(chr, chr1) == 0 && pos < pos1){
          error("input sequence mapping resutls are not sorted\n");
        }
        
        if(k5 > 0 || k3 > 0){
          sprintf(line, "%s\t%ld\t%d\t%d\n", chr1, pos1, k5, k3);
          fputs(line, FO);
        }

        strcpy(chr1, chr);
        pos1 = pos;
        k5 = k3 = 0;
  
        if(strcmp(strand, "+")==0){
          k5++;
        }else if(strcmp(strand, "-")==0){
          k3++;
        }else{
          error("unrecognized strand %s \n", strand);
        }

      }   
    }
        
  }
  
  fclose(FI);
  fclose(FO);
}


