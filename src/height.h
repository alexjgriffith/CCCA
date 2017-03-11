/* This file is part of peakAnalysis,
   http://github.com/alexjgriffith/alpha-score/, 
   and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
   the three-clause BSD License; see LICENSE.txt.
   Author : Alexander Griffith
   Contact: griffitaj@gmail.com */

//#ifndef CCCA_
//#include "CCCA.h"
//#define CCCA_
//#endif 
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>



#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int getChromosome(char value[10]);
int getChromosomeShort(char value[2]);
void getChromosomeValue(int order,char * value);
void rankChromosomes(char ** chroms,int *length, int *out);
void valueChromosomes(char ** chroms,int *length, int *out);
int compare(int start1, int end1,int start2, int end2);
void pileup(char ** filename,char ** chro,int *start,
	    int *end,int *length,int *scores);
void file_length(char ** filename,int * i);
void read_bed(char ** filename,char ** chrom,int *start, int *end);
void getChroms(char ** filename,char ** chroms);


typedef struct chromosomes {
  int order;
  char value[10];
  char shortValue[2];} chromosomes;

// The chromosome order output of bwa is unusual
static const chromosomes hg19Chrom[] = {
{0,"chr1","1"},
{1,"chr2","2"},
{2,"chr3","3"},
{3,"chr4","4"},
{4,"chr5","5"},
{5,"chr6","6"},
{6,"chr7","7"},
{7,"chrX","X"},
{8,"chr8","8"},
{9,"chr9","9"},
{10,"chr10","10"},
{11,"chr11","11"},
{12,"chr12","12"},
{13,"chr13","13"},
{14,"chr14","14"},
{15,"chr15","15"},
{16,"chr16","16"},
{17,"chr17","17"},
{18,"chr18","18"},
{19,"chr20","20"},
{20,"chrY","Y"},
{21,"chr19","19"},
{22,"chr22","22"},
{23,"chr21","21"},
{24,"chrM","M"},
{25,"chr",""}};

