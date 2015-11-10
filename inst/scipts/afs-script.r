#!/usr/bin/env R
#
# This file is part of CCCA,
# http://github.com/alexjgriffith/CCCA/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com

library(CCCA)
library(parallel)

# categories<-c("cd34_new","eryt","jurk","cem_1")
# categories<-as.character(unlist(read.table("~/Dropbox/UTX-Alex/jan/catagories")))
# peakXLSFiles<-paste("~/Dropbox/Data/august_peaks/",categories,"~combined_mock_peaks.xls",sep="")
#a<-makeAFS(peakXLSFiles,categories)
# rawReadFiles<-paste(categories,"/",categories,"_unique_nodupes.bed",sep="")
# print(dim(a))
#peakBedFiles<-paste("~/.peaktemp/",categories,"~combined",sep="")
#Rprof("inst/rprof/test-bed.r")
#a<-makeAFS(peakBedFiles,categories,"bed") 
#Rprof(NULL)

setwd("/mnt/brand01-00/mbrand_analysis")
    
for(pvalue in c(0,5,7.5,10,12.5,15)){
    categories<-c("cd34_new","eryt","jurk","cem_1")
    peakXLSFiles<-paste("peaks/october/",categories,"/combined_mock_peaks.xls",sep="")
    rawDataFiles<-paste("data_sets/",categories,"/",categories,"_unique_nodupes.bed",sep="")
    a<-makeAFS(peakXLSFiles,categories,pValue=pvalue)
    width<-dim(a)[2]
    colnames(a)<-c("chro",colnames(a)[2:width])
    a<-hg19Sort(a)
    a<-a[a$chro!="chrY"& a$chro!="chrM"& a$chro!="chrY",]
    score<-pileUp(a,rawDataFiles,4,TRUE)
    colnames(score)<-categories

    frontName<-function(x) paste("~/thesis-november/4x4-",pvalue,"-",x,".data",sep="")
    write.table(1-cor(score),frontName("cor"),quote=FALSE,sep="\t")
    write.table(1-cor(qn(score)),frontName("qn-cor"),quote=FALSE,sep="\t")
    write.table(1-overlap(a[,4:width]),frontName("ovr"),quote=FALSE,sep="\t")
    write.table(a,frontName("score-matrix"),quote=FALSE,row.names=FALSE,sep="\t")
    write.table(score,frontName("overlap-matrix"),quote=FALSE,row.names=FALSE,sep="\t")
    
    categories<-c("k562_1","eryt","jurk","cem_1")
    peakXLSFiles<-paste("peaks/october/",categories,"/combined_mock_peaks.xls",sep="")
    rawDataFiles<-paste("data_sets/",categories,"/",categories,"_unique_nodupes.bed",sep="")
    a<-makeAFS(peakXLSFiles,categories,pValue=pvalue)
    width<-dim(a)[2]
    colnames(a)<-c("chro",colnames(a)[2:width])
    a<-hg19Sort(a)
    a<-a[a$chro!="chrY"& a$chro!="chrM"& a$chro!="chrY",]
    score<-pileUp(a,rawDataFiles,4,TRUE)
    colnames(score)<-categories
    frontName<-function(x) paste("~/thesis-november/4x4-b-",pvalue,"-",x,".data",sep="")
    write.table(1-cor(score),frontName("cor"),quote=FALSE,sep="\t")
    write.table(1-cor(qn(score)),frontName("qn-cor"),quote=FALSE,sep="\t")
    write.table(1-overlap(a[,4:width]),frontName("ovr"),quote=FALSE,sep="\t")
    write.table(a,frontName("score-matrix"),quote=FALSE,row.names=FALSE,sep="\t")
    write.table(score,frontName("overlap-matrix"),quote=FALSE,row.names=FALSE,sep="\t")


    overlap<-abind(lapply(c(0,5,7.5,10,12.5,15),function(x)read.table(paste("4x4-",x,"-ovr.data",sep=""))),along=3)
    correlation<-abind(lapply(c(0,5,7.5,10,12.5,15),function(x)read.table(paste("4x4-",x,"-cor.data",sep=""))),along=3)

    plot(correlation[3,4,]-correlation[1,2,],type="l")
    
    plot(overlap[3,4,]-overlap[1,2,],type="l")
}
