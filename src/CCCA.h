/* This file is part of CCCA,
   http://github.com/alexjgriffith/CCCA/, 
   and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
   the three-clause BSD License; see LICENSE.txt.
   Author : Alexander Griffith
   Contact: griffitaj@gmail.com */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#define R_NO_REMAP
#define MAX_BED_BUFFER_SIZE 1024

#include <R.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Error.h>
#include <Rinternals.h>

void pileup(char ** filename,char ** chro,int *start,
	    int *end,int *length,int *scores);

void file_length(char ** filename,int * i);

void read_bed(char ** filename,char ** chrom,int *start, int *end);

void unityOutput(int * intChr, int * intSummit, int * name, int * l1,
		 int * l2, int * peaklength, int* peakwidth,
		 int* retChr, int * retSummit, int * retMatrix);


static R_NativePrimitiveArgType read_bed_t[]={STRSXP,STRSXP,INTSXP,INTSXP};
static R_NativePrimitiveArgType file_length_t[]={STRSXP,INTSXP};
static R_NativePrimitiveArgType pileup_t[]={STRSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP};
static R_NativePrimitiveArgType unityOutput_t[]={INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP};


static R_CMethodDef cMethods[]={
  {"read_bed",(DL_FUNC) &read_bed,4,read_bed_t},
  {"file_length",(DL_FUNC) &file_length,2, file_length_t},
  {"pileup",(DL_FUNC) &pileup,6, pileup_t},
  {"unityOutput",(DL_FUNC) &unityOutput,10,unityOutput_t},
  {NULL,NULL,0}
};

void R_init_ccca(DllInfo *info);
