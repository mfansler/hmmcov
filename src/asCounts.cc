/**********************************************************************
 *
 * asCounts.c
 *
 * allele specfic sequecing study
 *
 * copyright (c) 2009, Wei Sun, UNC-CH
 *
 * last modified May 04, 2009
 * first written May 04, 2009
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 **********************************************************************/
#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include "stdhash.hh"

#define MAX_LEN 1023
#define ZERO_QU 33

extern "C" {
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

hash_map_char<char*> *load_snp_set(FILE *fp)
{
	char chr[255], allele1[7], allele2[7];
	char buffer[2], key[255];
	int pos, c;
	hash_map_char<char*> *hash = new hash_map_char<char*>;
	while (fscanf(fp, "%s\t%d\t%s\t%s", chr, &pos, allele1, allele2) == 4) {
		sprintf(key, "%s.%d", chr, pos);
		buffer[0] = allele1[0];
		buffer[1] = allele2[0];
		buffer[2] = 0;
		while ((c = fgetc(fp)) != EOF && c != '\n');
		hash->insert(key, strdup(buffer));
	}
	return hash;
}

void asCounts(char **Rinput, char **Routput, char ** RsnpList, 
 int *Rm, int *RQ, int *Rq, int* RsnpQ, int *Rnl0, int *Rnu0, 
 int *Rn1, int *skip)
{
/**
 *
  m INT 	Maximum number of mismatches allowed for a read to be used [4]
  Q INT 	Maximum allowed number of quality values of mismatches [40]
  q INT 	Minimum mapping quality allowed for a read to be used [1] 
  nl0 INT Maximum allowed number of 0-mismatch hits of the first 24bp [0]
  nu0 INT Maximum allowed number of 0-mismatch hits of the first 24bp [85]
  n1  INT Maximum allowed number of 1-mismatch hits of the first 24bp [85]
 */
  int m = *Rm;
  int Q = *RQ;
  int q = *Rq;
  int snpQ = *RsnpQ;
  int nl0 = *Rnl0;
  int nu0 = *Rnu0;
  int n1  = *Rn1;
  int i, j, maq, nm, sqm, k50, k51, k30, k31, nm0, nm1, L;
  long int pos, pos1;
  FILE *FI, *FO, *Fsnp;
  char str_buf[MAX_LEN], line[MAX_LEN];
  char *key, *alleles;
  char chr[127], chr1[127], strand[2], seq[255], qua[255];
	const char *delim = "\t";
	char *ptr= NULL;
  char *input, *output, *snpList;
  input   = Rinput[0];
  output  = Routput[0];
  snpList = RsnpList[0];
  
  int n_unrecognized = 0;
  
	key = (char*)calloc(255, 1);
	alleles = (char*)calloc(7, 1);

  // Rprintf("snpQ + ZERO_QU=%d\n", snpQ+ZERO_QU);

  /**
   * read in the snp information
   */
  Fsnp = fopen(snpList, "r");
  assert(Fsnp);
	hash_map_char<char*> *hash_map = 0;
  hash_map = load_snp_set(Fsnp);
  fclose(Fsnp);

  // Rprintf("input is %s, and output is %s\n", input, output);
  
  FI = fopen(input, "r");
  assert(FI);
  FO = fopen(output, "w");
  assert(FO);
  
  for(i=0; i < *skip; i++){
    if(fgets(str_buf, MAX_LEN, FI) == NULL){
      error("thre are only %d lines in file %s\n", i, input);
    }
  }

  /**
   * read in the first line
   * obtain the initial chromosome (chr1) and position (pos1)
   */
  k50 = k51 = k30 = k31 = 0;
  
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
      }else if(i==14){
        L   = atoi(ptr);
      }else if(i==15){
        strcpy(seq, ptr);
      }else if(i==16){
        strcpy(qua, ptr);
      }
    } while ((ptr = strtok(NULL, delim)) != NULL);
  }else{
    error("%s is not tab-delimated\n", input);
  }

  // Rprintf("chr: %s, pos: %d, strand: %s, length: %d, ", chr1, pos1, strand, L);
  // Rprintf("score: %d, nmismatch: %d, score of mismatch: %d\n", maq, nm, sqm);
  
  if(maq >= q && nm <= m && sqm <= Q && nm0 <= nu0 && nm0 >= nl0 && nm1 <= n1){
    for(j=0; j<L; j++){
  		/**
  		 * make sure the quality of this particular base is no less than snpQ
  		 */
  		if((int)qua[j] < ZERO_QU + snpQ){
  		  continue;
  		}
  		
  		sprintf(key, "%s.%ld", chr, pos+j);
      if (hash_map->find(key, &alleles)) {
                
        if(seq[j] == alleles[0]){
          if(strcmp(strand, "+")==0){
            k50++;
          }else if(strcmp(strand, "-")==0){
            k30++;
          }else{
            error("unrecognized strand %s \n", strand);
          }
        }else if(seq[j] == alleles[1]){
          if(strcmp(strand, "+")==0){
            k51++;
          }else if(strcmp(strand, "-")==0){
            k31++;
          }else{
            error("unrecognized strand %s \n", strand);
          }
        }
        break;
      }
    }
  }

  // Rprintf("read in other lines\n");

  /**
   * read in the other lines
   */

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
        }else if(i==14){
          L   = atoi(ptr);
        }else if(i==15){
          strcpy(seq, ptr);
        }else if(i==16){
          strcpy(qua, ptr);
        }
      } while ((ptr = strtok(NULL, delim)) != NULL);
    }else{
      error("%s is not tab-delimated\n", input);
    }
    
    // Rprintf("chr: %s, pos: %d, strand: %s, length: %d, ", chr, pos, strand, L);
    // Rprintf("seq: %s\n", seq);
    // Rprintf("qua: %s\n", qua);
    
    if(maq >= q && nm <= m && sqm <= Q && nm0 <= nu0 && nm0 >= nl0 && nm1 <= n1)
    {
      if(strcmp(chr, chr1) == 0 && pos == pos1){
                
        for(j=0; j<L; j++){
          
          /**
           * make sure the quality of this particular base is no less than snpQ
           */
          if((int)qua[j] < ZERO_QU + snpQ){
            continue;
          }
          
          sprintf(key, "%s.%ld", chr, pos+j);
          if (hash_map->find(key, &alleles)) {
          
            if(seq[j] == alleles[0]){
              if(strcmp(strand, "+")==0){
                k50++;
              }else if(strcmp(strand, "-")==0){
                k30++;
              }else{
                error("unrecognized strand %s \n", strand);
              }
            }else if(seq[j] == alleles[1]){
              if(strcmp(strand, "+")==0){
                k51++;
              }else if(strcmp(strand, "-")==0){
                k31++;
              }else{
                error("unrecognized strand %s \n", strand);
              }
            }
            break;
          }
        }
      }else{
        if(strcmp(chr, chr1) == 0 && pos < pos1){
          error("input sequence mapping resutls are not sorted\n");
        }
        
        if(k50 > 0 || k30 > 0 || k51 > 0 || k31 > 0){
          sprintf(line, "%s\t%ld\t%d\t%d\t%d\t%d\n", chr1, pos1, k50, k30, k51, k31);
          fputs(line, FO);
        }
                
        strcpy(chr1, chr);
        pos1 = pos;
        k50 = k51 = k30 = k31 = 0;
                
        for(j=0; j<L; j++){
          /**
           * make sure the quality of this particular base is no less than snpQ
           */
          if((int)qua[j] < ZERO_QU + snpQ){
            continue;
          }
                    
          sprintf(key, "%s.%ld", chr, pos+j);
          if(hash_map->find(key, &alleles)) {
                        
            if(toupper(seq[j]) == alleles[0]){
              if(strcmp(strand, "+")==0){
                k50++;
              }else if(strcmp(strand, "-")==0){
                k30++;
              }else{
                error("unrecognized strand %s \n", strand);
              }
            }else if(toupper(seq[j]) == alleles[1]){
              if(strcmp(strand, "+")==0){
                k51++;
              }else if(strcmp(strand, "-")==0){
                k31++;
              }else{
                error("unrecognized strand %s \n", strand);
              }
            }else {
              n_unrecognized += 1;
              warning("observed allele %c, expect %c or %c\n", seq[j], alleles[0], alleles[1]);
            }

            break;
          }
        }

      }
    }
 
  }
  
  Rprintf("there are %d un-recognized SNP locations\n", n_unrecognized);
  
  fclose(FI);
  fclose(FO);
  delete hash_map;
	free(key);
	free(alleles);
}
}// extern "C"
