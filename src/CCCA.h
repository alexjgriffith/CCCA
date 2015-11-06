#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

// heights.c
int getChromosome(char value[10]);
int getChromosomeShort(char value[2]);
void getChromosomeValue(int order,char * value);
void rankChromosomes(char ** chroms,int *length, int *out);
void valueChromosomes(char ** chroms,int *length, int *out);
int compare(int start1, int end1,int start2, int end2);
void pileup(char ** filename,int * chro,int *start,int *end,int *peaknum,int *scores);
void file_length(char ** filename,int * i);
void read_bed(char ** filename,char ** chrom,int *start, int *end);
void getChroms(char ** filename,char ** chroms);

// regions.c
int addToBuffer(int logical, int category, int *buffer);
void region(int * value, int * logical , int * category,int * buffer, int * outMatrix,int * length,int * width);
void unityOutput(int * intChr, int * intSummit, int * name, int * l1, int * l2, int * peaklength, int* peakwidth, int* retChr, int * retSummit, int * retMatrix);

#ifndef CME_
static R_NativePrimitiveArgType read_bed_t[]={STRSXP,STRSXP,INTSXP,INTSXP};
static R_NativePrimitiveArgType file_length_t[]={STRSXP,INTSXP};
static R_NativePrimitiveArgType pileup_t[]={STRSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP};
static R_NativePrimitiveArgType rankChromosomes_t[]={STRSXP,INTSXP,INTSXP};
static R_NativePrimitiveArgType valueChromosomes_t[]={STRSXP,INTSXP,INTSXP};
static R_NativePrimitiveArgType getChroms_t[]={STRSXP,STRSXP};
static R_NativePrimitiveArgType region_t[]={INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP};

static R_NativePrimitiveArgType unityOutput_t[]={INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP,INTSXP};


static R_CMethodDef cMethods[]={
{"read_bed",(DL_FUNC) &read_bed,4,read_bed_t},
{"file_length",(DL_FUNC) &file_length,2, file_length_t},
{"pileup",(DL_FUNC) &pileup,6, pileup_t},
{"rankChromosomes",(DL_FUNC) &rankChromosomes,6,rankChromosomes_t},
{"valueChromosomes",(DL_FUNC) &valueChromosomes,6,valueChromosomes_t},
{"getChroms",(DL_FUNC) &getChroms,3,getChroms_t},
{"region",(DL_FUNC) &region,7,region_t},
{"unityOutput",(DL_FUNC) &unityOutput,10,unityOutput_t},
{NULL,NULL,0}
  };
#define CME_
#endif

 void R_init_myLib(DllInfo *info);
