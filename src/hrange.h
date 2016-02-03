#include <stdio.h>
#include <stdlib.h>

#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>



typedef struct peak{
  char chr[100];
  int chr_ord;
  int start;
  int end;
} peak;


void rheight(char * filename,peak * peaks,int length, peak *** scores,int ** heights);

int convertHeights(peak * temp, int  length, peak ** scores, int * lengths,int *** collectIn);

void  peakDensity(char ** filename,char ** chro,int *start,
		  int *end,int *length,int *scoresOut);

void buildPeaks(char **chr, int * start, int * end , int length, peak ** temp);
