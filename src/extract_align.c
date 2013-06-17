/**********************************************************************
 *
 * alignbility
 *
 * extract alignbility at specificied locations
 * 
 *
 * copyright (c) 2009, Wei Sun, UNC-CH
 *
 * last modified Dec 04, 2009
 * first written Dec 04, 2009
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
 * extract alignbility at specificied locations
 * For each location (bp), extract alignbility for a region of
 * of certain length (sequence tag length) from one side
 * (depens on strand)
 *
 * There are two input files:
 *
 * One file includes the base pair locations where there are
 * at least one base count:
 *
 * track type=wiggle_0 name="aa" desc="bb" visibility=full
 * variableStep chrom=chrX
 * 3824282 1
 * 3824283 1
 * 3824284 1
 * 3824285 1
 * 
 * The other file includes the genome-wide alignbility:
 * 
 * >chrX
 * 00000000000000000000000000000000000000000000000000
 * 00000000000000000000000000000000000000000000000000
 *
 **********************************************************************/

void extract_align(char **Rinput1, char **Rinput2, char **Routput, int* dim)
{

  char *input1, *input2, *output;
  int skip1, skip2, i, j, pos1, count1, nrow1;
  int align1, tagLen, readLen, strand, ct1, chrLen;
  char str_buf[MAX_LEN], line[MAX_LEN];
  long int * posV;
  long int start;
  long int end;
  short int * alignV;
  FILE *FI1, *FI2, *FO;

  input1 = Rinput1[0];
  input2 = Rinput2[0];
  output = Routput[0];
  
  nrow1  = dim[0];
  skip1  = dim[1];
  skip2  = dim[2];
  chrLen = dim[3];
  tagLen  = dim[4] - 1;
  readLen = dim[5] - 1;
  strand  = dim[6];
  
  posV = (long int *)Calloc(nrow1, long int);
  alignV = (short int *)Calloc(chrLen+1, short int);

  Rprintf("input1: %s\ninput2: %s\noutput: %s\n", input1, input2, output);
  
  FI1 = fopen(input1, "r");
  if(FI1 == NULL){
    error("cannot open file %s\n", input1);
  }
  
  FI2 = fopen(input2, "r");
  if(FI2 == NULL){
    error("cannot open file %s\n", input2);
  }
    
  FO = fopen(output, "w");
  if(FO == NULL){
    error("cannot open file %s\n", output);
  }
  
  /**
   * skip first a few lines for both files
   */
  
  for(i=0; i < skip1; i++){
    if(fgets(str_buf, MAX_LEN, FI1) == NULL){
      error("thre are only %d lines in file %s\n", i, input1);
    }
  }

  for(i=0; i < skip2; i++){
    if(fgets(str_buf, MAX_LEN, FI2) == NULL){
      error("thre are only %d lines in file %s\n", i, input2);
    }
  }
  
  /**
   * read the first file, obtain the locations
   */
  
  i = 0;
  while(fscanf(FI1, "%d %d", &pos1, &count1) != EOF){
    posV[i] = pos1;
    i++;
  }
    
  Rprintf("finish reading %s, i=%d\n", input1, i);
  if (i != nrow1) {
    error("file %s has %d rows, not %d rows\n", input1, i, nrow1);
  }
  
/*
  for (i=0; i<100; i++) {
    Rprintf("%d\n", posV[i]);
  }
*/
  
  /**
   * read the second file, obtain the alignbility
   */
  i = 0;
  align1 = fgetc(FI2);
  while( align1 != EOF){
    align1 = (short int) align1 - 48;
    
    if (align1 >= 0) {
      i++;
      alignV[i] = align1;
    }
    
    align1 = fgetc(FI2);
  }

  Rprintf("finish reading %s, i=%d\n", input2, i);
/*
  for (i=0; i<100; i++) {
    Rprintf("%d\n", alignV[i]);
  }
*/
  fclose(FI1);
  fclose(FI2);
  
  /**
   * for each position in posV, calculate its alignbility
   */
  
  for (i=0; i<nrow1; i++) {
    pos1 = posV[i];
    ct1  = 0;

    if (strand >= 0) {
      start = pos1 - tagLen;
      end   = pos1;
      if(start <= 1){ start = 1; }
      
      for (j=start; j<=end; j++) { ct1 += alignV[j]; }
    }
    
    if (strand <= 0) {
      start = pos1 - readLen;
      end   = pos1 + tagLen - readLen;
      if(start < 1){ start = 1; }
      if(end > chrLen){ end = chrLen; }
      
      for (j=start; j<=end; j++) { ct1 += alignV[j]; }
    }
    
    fprintf(FO, "%d\n", ct1);
  }
  
  fclose(FO);

  free(posV);
  free(alignV);
}


