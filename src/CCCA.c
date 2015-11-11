/* This file is part of CCCA,
   http://github.com/alexjgriffith/CCCA/, 
   and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
   the three-clause BSD License; see LICENSE.txt.
   Author : Alexander Griffith
   Contact: griffitaj@gmail.com */

#include "CCCA.h"

void R_init_myLib(DllInfo *info)
{
R_registerRoutines(info,cMethods,NULL,NULL,NULL);
}
