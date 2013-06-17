/**********************************************************************
 *
 * base counts
 *
 * count the number of overlapped reads at each base pair location
 * 
 *
 * copyright (c) 2009, Wei Sun, UNC-CH
 *
 * last modified Dec 05, 2009
 * first written Dec 05, 2009
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
 *
 * There is one input files:
 *
 * One file includes the sequence mapping results at one chromosome
 * with two columns, one is strand, and the other one is location
 *
 *  -       91009338
 *  -       100662078
 *  -       132749072
 *  -       66858796
 *  +       17855825
 *  +       12700392
 *  -       100673664
 *  +       39460698
 * 
 **********************************************************************/

void baseCounts(char **Rinput1, char **Routput, int* dim)
{
  char *input1, *output;
  int skip, i, *ctV, chrLen, readLen, tagLen, start, end;
  long int pos1;
  char str_buf[MAX_LEN], strand1[MAX_LEN];
  
  FILE *FI1, *FO;

  input1 = Rinput1[0];
  output = Routput[0];
  
  skip    = dim[0];
  
  // assign memory for one more number, since index start from 1
  chrLen  = dim[1] + 1;
  
  // subtract 1 to simplifiy the caculation since: end = start + Len - 1
  readLen = dim[2] - 1;
  tagLen  = dim[3] - 1;

  ctV = (int *)Calloc(chrLen, int);

  // Rprintf("input1: %s, output: %s\n", input1, output);
  // Rprintf("chrLen=%d, readLen=%d, tagLen=%d\n", chrLen, readLen, tagLen);

  FI1 = fopen(input1, "r");
  if(FI1 == NULL){
    error("cannot open file %s\n", input1);
  }
      
  FO = fopen(output, "w");
  if(FO == NULL){
    error("cannot open file %s\n", output);
  }
  
  /**
   * skip first a few lines for both files
   */
  
  for(i=0; i < skip; i++){
    if(fgets(str_buf, MAX_LEN, FI1) == NULL){
      error("thre are only %d lines in file %s\n", i, input1);
    }
  }
  
  /**
   * read the first file, obtain the strands and counts
   */
  while(fscanf(FI1, "%s %ld", strand1, &pos1) != EOF){
    if(strand1[0] == '+'){
      start = pos1;
      end   = start + tagLen;
      if (end > chrLen) { end = chrLen; }
      for (i=start; i<=end; i++) {
        ctV[i] += 1;
      }
    }else if (strand1[0] == '-') {
      start = pos1 + readLen;
      end   = start - tagLen;
      if (end < 1) { end = 1; }
      for (i=start; i>=end; i--) {
        ctV[i] += 1;
      }      
    }
  }
  
  for (i=1; i<=chrLen; i++) {
    if (ctV[i] > 0) {
      fprintf(FO, "%d %d\n", i, ctV[i]);
    }
  }
  
  fclose(FI1);
  fclose(FO);

  Free(ctV);
}


