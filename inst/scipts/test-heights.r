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
install()
library(CCCA)

#getChroms("~/Dropbox/UTX-Alex/jan/combined_mock.bed")

cats<-read.table("/home/griffita/Dropbox/UTX-Alex/jan/catagories")
prefix<-"/mnt/brand01-00/mbrand_analysis/data_sets/"
suffix<-"_sorted.bed"
rawdata<-apply(cats,1,function(x){paste(prefix,x,"/",x,suffix,sep="")})
# Apply pileUp to peaks found using MACS
data<-hg19Sort(loadBedFile(file))
score<-pileUp(data,rawdata,n=22)
 
