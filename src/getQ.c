/**********************************************************************
 *
 * getQ.c
 *
 * simply change the ASCII code to the numeric value
 *
 * copyright (c) 2010, Wei Sun, UNC-CH
 *
 * last modified Feb 27, 2010
 * first written Feb 27, 2010
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

#define MAX_LEN 1023

/**********************************************************************
 *
 * simply change the ASCII code to its numeric value
 *
 **********************************************************************/

void getQ(char **Rinput, int* output, int* n)
{
  int i;
  char *input;
  input   = Rinput[0];
  
  for(i=0; i<*n; i++){
    output[i] = (int) input[i];
  }
  
}
