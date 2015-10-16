#!/usr/bin/env R
#
# This file is part of CCCA,
# http://github.com/alexjgriffith/CCCA/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com


library(devtools)
setwd("~/Masters/CCCA")
document()
library(CCCA)
getChroms("~/Dropbox/UTX-Alex/jan/combined_mock.bed")
