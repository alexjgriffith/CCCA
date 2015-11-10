#!/usr/bin/env R
#
# This file is part of CCCA,
# http://github.com/alexjgriffith/CCCA/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com


#' unify bed file
#'
#' unify bed file
#' 
#' @export
unifyBedFile<-function(dataset,overlapWidth, chr="chr",summit="summit",name="name"){
    chrcats<-levels(dataset[chr,])
    l<-nrow(dataset)
    t1<-dataset[1:l-1,]
    t2<-dataset[2:l,]
    b<-(t1[,chr]==t2[,chr] & abs(as.numeric(t1[,summit])-as.numeric(t2[,summit]))<overlapWidth)
    prePeaks<-unlist(lapply(which(b==FALSE),function(x){c(x,x+1)}))
    peaks<-t(matrix(c(1,prePeaks[1:length(prePeaks)-1]),nrow=2))
    data<-unityOutput(peaks,
                as.integer(dataset[,chr]),
                as.integer(dataset[,summit]),
                as.integer(dataset[,name]))
    names(data)<-list(chr,summit,"matrix")
    data[[chr]]<-as.factor(data[[chr]])
    levels(data[[chr]])<-levels(dataset[,chr])
    colnames(data$matrix)<-levels(dataset[,name])
    od<-data.frame(chr=data[[chr]],summit=data[[summit]])
    od<-cbind(od,data$matrix)
    od
}

#' unity output
#'
#' Wrapper for the unityOutput C function.
#' 
#' @param peaks A 2xN matrix representing overlaping peak regions
#' @param intChr N length list of chomosomes represented as integers
#' @param intSummit The summit location 
#' @param intname The data set or name that the peak come froms in int form
#'  @export
unityOutput<-function(peaks,intChr,intSummit,intname){
    nchr<-length(unique(intname))
    lpeaks<-dim(peaks)[1]
    retChr<-integer(lpeaks)
    retMatrix<-integer(lpeaks*nchr)
    retSummit<-integer(lpeaks)
    data<-.C("unityOutput",as.integer(intChr),as.integer(intSummit),as.integer(intname),as.integer(peaks[,1]),as.integer(peaks[,2]),as.integer(lpeaks),as.integer(nchr),chr=retChr,summit=retSummit,matrix=retMatrix)
    list(chro=data$chr,summit=data$summit,matrix=t(matrix(data$matrix,nrow=nchr)))
   }


#' Shift from zero
#' 
#' Checks if the result of subtracting width from value is zero
#'
#' @param value integer or double which cannot be less than 0
#' @param width integer or double
#' @return a range value+width value-width where value-width >0
#' @template authorTemplate

#' @export
shiftFromZero<-function(summit){
   testBedSE<-cbind(summit-350,summit+350)
   x<-which(testBedSE[,1]<0)
   testBedSE[x,2]=testBedSE[x,2]-testBedSE[x,1]
   testBedSE[x,1]=0
   testBedSE
}


#' Make AFS
#'
#' 
#' @export
makeAFS<-function(peakList,categories,format="xls",pValue=FALSE){
    if(format!="xls" & pValue){
        print("to have a pvalue cut off the file format must be in MACS XLS output")
        return(0)
    }
    readPeaksXLS<-function(file,name=file){
        bedData<-read.table(file,header=TRUE,skip="#")
        bedData<-bedData[bedData$X.log10.pvalue>pValue,]
        ret<-data.frame(bedData$chr,bedData$abs_summit,name)
        colnames(ret)<-c("chr","summit","name")
        ret
    }
    readPeaksBed<-function(file,name=file){
        # much faster, but far more error prone
        bedData<-loadBedFile(file)
        colnames(bedData)<-c("chr","start","end")
        summit<-bedData$start+floor((bedData$end-bedData$start)/2)
        cbind(bedData,name,summit)
    }
    loadBedFun<-switch(format,
                       xls=readPeaksXLS,
                       bed=readPeaksBed)
    testBed<-unifyBedFile(
        sortDataFrame(do.call(rbind,lapply( mapziplist(peakList,categories),function(x) do.call(loadBedFun, as.list(x)))),"chr","summit"),700)
    return(testBed)
    shiftFromZero(testBed$summit)
    width<-dim(testBed)[2]
    retData<-data.frame(chr=testBed[,1], start=testBedSE[,1],end=testBedSE[,2])
    #retData<-data.frame(chr=testBed[,1], summit=testBed$summit)
    cbind(retData,testBed[,3:width])
       
}
