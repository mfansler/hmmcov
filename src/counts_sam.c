/**********************************************************************
 *
 * counts_sam.c
 *
 * allele specfic sequecing study
 *
 * copyright (c) 2010, Wei Sun, UNC-CH
 *
 * first written Aug 22, 2010
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
 * generate counts from the ouptut of tophat: accepted_hits.sam
 *
 **********************************************************************/

/**
 *

 Three lines of the input file looks like the following:
 
 IL26_1382:1:109:660:632 137     chr1    309     1       36M     *       0       0
    CAACCCCAACCCTAACCCTAACCCCTAACCCTAACC    <<=<<<<<<<4<<<9<;<<9<<;<;<;7<;<9:0:<
    NM:i:1  NH:i:4  CC:Z:chr2       CP:i:114076909
 
 IL26_1382:1:61:717:749  99      chr1    3037    3       36M     =       3093    0
    TGGGGAGGCAGCTGTAACTCAAAGCCTTAGCCTCTG    >>>>><>>>3>>;>:96=>9>2>><:8<39:+1::5
    NM:i:0  NH:i:2  CC:Z:chr9       CP:i:3287
 
 IL26_1382:1:61:717:749  147     chr1    3093    3       36M     =       3037    0
    TCAGGCGCCAAAGGGATTCTGCCAGCATAGTGCTCC    3;;&<6.><<<>><>7<<<<><>>>><<<>><<>>>
    NM:i:1  NH:i:2  CC:Z:chr9       CP:i:3343

 Each line of maq mapview output consists of
 [1] QNAME
 [2] FLAG
 [3] RNAME   reference name (chromsome)
 [4] POS     leftmost position of the sequence
 [5] MAPQ    phred-scaled quality score
 [6] CIGAR   M: match
 [7] MRNM    mate reference sequence
 [8] MPOS    left most position of mate sequence
 [9] ISIZE   inferred insert size
 [10] SEQ    Sequence
 [11] QUAL   ASCII - 33 gives Phred base quality
 [12] TAG:VTYPE:VALUE
             (1) NM: number of nucleotide difference
             (2) NH: number of alignments
             (3) CC: reference name of next hit
             (4) CP: Leftmost coordinate of the next hit
 */

void counts_sam(char **Rinput, char **Routput, int *Rflag2kp, int *Rflag2rm,  
 int *Rmin_mapQ, int *skip, int *simplify)
{
/**
 *
  flag2kp  INT 	flags to keep
  flag2rm  INT 	flags to remove
  min_mapQ INT 	Minimum mapping quality allowed for a read to be used [0] 
 */
  
  int flag2kp  = *Rflag2kp;
  int flag2rm  = *Rflag2rm;
  int min_mapQ = *Rmin_mapQ;
  
  long int pos1, pos2;
  
  int i, flag, mapQ=0, keepIt=1;
  
  FILE *FI, *FO;
  
  char str_buf[MAX_LEN], line[MAX_LEN];
  char chr[128];
  
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

    
  while(fgets(str_buf, MAX_LEN, FI) != NULL){
    
    keepIt = 1;
    
    if(str_buf[0] == '@') continue;
    
    strcpy(line, str_buf);
    
    i = 0;
    
    if ((ptr = strtok(str_buf, delim)) != NULL) {
      
      do {
        i++;
        if(i==2){
          flag = atol(ptr);

          if ((flag & flag2kp != flag2kp) || (flag & flag2rm)){
            keepIt = 0;
            break;
          }

        }else if(i==3){
          strcpy(chr, ptr);
        }else if(i==4){
          pos1 = atol(ptr);
        }else if(i==5){
          mapQ = atoi(ptr) - 33;
        }else if(i==8){
          pos2 = atol(ptr);
        }
      } while ((ptr = strtok(NULL, delim)) != NULL);
      
    }else{
      error("%s is not tab-delimated\n", input);
    }
    
    if(keepIt && mapQ >= min_mapQ)
    {
      if (*simplify) {
        sprintf(line, "%s\t%ld\t%ld\n", chr, pos1, pos2);
        fputs(line, FO);
      }else {
        fputs(line, FO);
      }
      
    }
  }
  
  fclose(FI);
  fclose(FO);
}
