/**********************************************************************
 *
 * map2SNP.cc
 *
 * Given a set of windows of fixed length, say 25bp
 * find the subset that harbor at least one SNP
 *
 * copyright (c) 2010, Wei Sun, UNC-CH
 *
 * last modified Feb 26, 2010
 * first written May 26, 2010
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
#include "stdhash.hh"

#define MAX_LEN 1023
#define ZERO_QU 33

extern "C" {

  hash_map_char<char*> *load_snp_loc(int *posSNP, int n)
  {
    int i;
    char key[255];
    hash_map_char<char*> *hash = new hash_map_char<char*>;

    for (i=0; i<n; i++) { 
      sprintf(key, "%d", posSNP[i]);
      hash->insert(key, strdup(key));
    }
    
    return hash;
  }

  void map2SNP(int* posSNP, int* posQuery, int* dims, int* found)
  {
    int i, j, pos;
    char *key, *value;
    key   = (char*)calloc(255, 1);
    value = (char*)calloc(255, 1);

    int nSNP   = dims[0];
    int nQuery = dims[1];
    int qLen   = dims[2];
    int nFound = 0;
    
    hash_map_char<char*> *hash_map = 0;
    hash_map = load_snp_loc(posSNP, nSNP);

    for (i=0; i<nQuery; i++) {
      //if(i % 1000 == 0){
      //  Rprintf("%d\n", i);
      //}
      
      pos = posQuery[i];
      nFound = 0;
      
      for(j=1; j<=qLen; j++){
        sprintf(key, "%d", pos+j);
        if (hash_map->find(key, &value)) { nFound += 1; }
      }
      
      found[i] = nFound;
    }
    
    delete hash_map;
    free(key);
    free(value);
    
  }
  
}// extern "C"
