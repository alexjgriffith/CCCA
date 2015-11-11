/* This file is part of CCCA,
   http://github.com/alexjgriffith/CCCA/, 
   and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
   the three-clause BSD License; see LICENSE.txt.
   Author : Alexander Griffith
   Contact: griffitaj@gmail.com */

#include <stdio.h>

#define R_NO_REMAP

#include <R.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Error.h>
#include <Rinternals.h>



int addToBuffer(int logical, int category, int *buffer);
void region(int * value, int * logical , int * category,int * buffer, int * outMatrix,int * length,int * width);
void unityOutput(int * intChr, int * intSummit, int * name, int * l1, int * l2, int * peaklength, int* peakwidth, int* retChr, int * retSummit, int * retMatrix);








