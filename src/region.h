//#ifndef CCCA_
//#include "CCCA.h"
//#define CCCA_
//#endif 


#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#include <stdlib.h>
#include <stdio.h>

int addToBuffer(int logical, int category, int *buffer);
void region(int * value, int * logical , int * category,int * buffer, int * outMatrix,int * length,int * width);
void unityOutput(int * intChr, int * intSummit, int * name, int * l1, int * l2, int * peaklength, int* peakwidth, int* retChr, int * retSummit, int * retMatrix);








