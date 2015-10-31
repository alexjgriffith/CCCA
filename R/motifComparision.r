#!/usr/bin/env R
#
# This file is part of CCCA,
# http://github.com/alexjgriffith/CCCA/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com
#
# functions:
# numberIntersect
# getHeights
# motifHist
# nearSummit
# getDistance
# motifs2View
# histVisualize
# heightHist
# locHist
# locHist2


### everything here needs a name overhaul

### Stand Alones

#' Intersect Count
#'
#' Determine the number of cases which 2 motifs intersect
#'
#' @export
numberIntersect<-function(motifA,motifB,Sequences,reg=rep(TRUE,length(Sequences))){
    #numberIntersect<-oldname
    motifALocations <- grep(unlist(IUPACtoBase(motifA)),Sequences)
    motifBLocations <- grep(unlist(IUPACtoBase(motifB)),Sequences)
    regionLocations <- which(reg)    
    length(intersect(intersect(motifALocations,motifBLocations)
                    ,regionLocations))
}

#' Get Heights
#'
#' Get the heights of the motif comparison data
#'
#' @export
getHeights<-function(h,range=c(min(h),max(h))){
    rep<-rep(0,(range[2]-range[1]+1))
    for(i in h)
        rep[i-range[1]+1]<-rep[i-range[1]+1]+1
    rep
}

#### Efficient comparision of motif locations

#' Motif histogram
#'
#' @examples
#' 
#' values<-loadHeightFile(heightFile)
#' data<-values$data
#' reg<-ascore(data,1,"top",3)
#' test<-readDNAStringSet(fastaFile,use.names = TRUE)
#' motifs<-as.matrix(read.table("data/normal_not_abnormal_motifs"))
#' mList<-unlist(lapply(c(motifs,addmotifs),IUPACtoBase))
#' cList<-unlist(lapply(lapply(c(motifs,"CGNNGC"),IUPACtoBase),compliment))
#' locationsM<-lapply(mList,grep,test)
#' locationsC<-lapply(cList,grep,test)
#' l<-length(cList)
#' motifHist(mList,cList,locationsM,locationsC,4,l,reg)
#' @export
 motifHist<-function(data,mList,cList,locationsM,locationsC,n1,n2,reg,width=c(-30,30)){
    lM<-intersect(intersect(locationsM[[n1]],locationsM[[n2]]),which(reg))
    lC<-intersect(intersect(locationsC[[n1]],locationsC[[n2]]),which(reg))
    bM<-c()
    bC<-c()
    if(length(lM)>0)
        bM<-lapply(c(mList[n1],mList[n2]),function(x) lapply(gregexpr(x, data[lM]),as.numeric))
    if(length(lC)>0)        
        bC<-lapply(c(cList[n1],cList[n2]),function(x) lapply(gregexpr(x, data[lC]),as.numeric))
    #h<-c(getDistance(bM[[1]],bM[[2]]),-getDistance(bC[[1]],bC[[2]]))
    #print(c(mList[n1],consensusIUPAC(mList[n2]))
    #print(nchar(mList[n1])-nchar(mList[n2]))
    print(nchar(consenusIUPAC(mList[n1]))-nchar(consenusIUPAC(mList[n2])))
    h<-c(getDistance(bM[[1]],bM[[2]]),-(nchar(consenusIUPAC(mList[n1]))-nchar(consenusIUPAC(mList[n2]))+getDistance(bC[[1]],bC[[2]])))
    #h<-getDistance(bM[[1]],bM[[2]])#getDistance(bM[[1]],bM[[2]])#,-getDistance(bC[[1]],bC[[2]]))
    #if(length(h)>0)
       # hist(h,breaks=1000,xlim=width,xlab=paste(mList[n1],mList[n2],sep=" - "))
    h}


#'@export
nearSummit<-function(data,mList,cList,locationsM,locationsC,n1,reg,width=150)
    {
    lM<-intersect(locationsM[[n1]],which(reg))
    lC<-intersect(locationsC[[n1]],which(reg))
    bM<-c()
    bC<-c()
    if(length(lM)>0){
        bM<-lapply(gregexpr(mList[n1], data[lM]),as.numeric)
        #print(mList[n1])
    }
    if(length(lC)>0)        {
        bC<-lapply(gregexpr(cList[n1], data[lC]),as.numeric)
        #print("b")
    }
    h<-c(getDistance(unlist(bM),width),getDistance(unlist(bC),width))
    h}

#' @export
getDistance<-function(x,y,one=FALSE){
    as.numeric(
        unlist(mapply(function(x,y){
        temp<-outer(x, y,"-")
        temp<-temp[upper.tri(x=temp,diag=TRUE)]
        #print(str(temp))
        if(one ){temp[which.max(abs(temp))]}
        else  {temp}
    },x,y)))}


### ineficent, but low overhead for comparion just a few data sets
#' Motifs Two View
#' 
#' Visualize the similarity between two motifs,
#' Very ineficent,
#' @export
#' @examples
#' FastaFile<-"~/Dropbox/UTX-Alex/jan/combined.fasta"
#' Sequences <- readDNAStringSet(FastaFile, "fasta")
#' data<-loadHeightFile("~/Dropbox/UTX-Alex/jan/combined_heights.bed")$data
#' reg<-as.matrix(ascore(data,c(1,1,1),c("top","bottom","middle"),c(1,1,1)))
#' m<-motifs2View("CANNTG","TGACCT",reg[,3],Sequences)
#' k<-motifs2View("CANNTG","GATAAG",reg[,1],Sequences)
#' t0<-unlist(lapply(seq(min(k),max(k)),function(i)length(which(k==i))))
motifs2View<-function(m1,m2,reg,Sequences,nearHeights=FALSE){
    m12<-c(m1, m2)
    mList<-unlist(lapply(m12,IUPACtoBase))
    # Reverse compliment
    cList<-lapply(unlist(lapply(m12,IUPACtoBase)),compliment)
    # Find locations of forward and reverse compilent
    locationsM<-lapply(mList,grep,Sequences)
    locationsC<-lapply(cList,grep,Sequences)
    if(nearHeights==FALSE){
        h<-list(motifHist(Sequences,
                          mList,
                          cList,
                          locationsM,
                          locationsC,1,2,reg))
        histVisualize(h,m1,m2)
        h<-unlist(h)
        t0<-unlist(lapply(seq(min(h),max(h)),function(i)length(which(h==i))))
    }
    else{
        h<-nearSummit(Sequences,mList,cList,locationsM,locationsC,1,reg)}
    #c((max(t0)-mean(t0))/sqrt(var(t0)),scoreFunction(t0))
    h
}


### Plotting Functions

#' @title Visualize batches of height data
#' @export
histVisualize<-function(h,m1,m2,n=1,name=c("NoName")){
    lengthN<-length(n)
    pushViewport(viewport(layout=grid.layout(lengthN,2)))    
    if(n==1){
        h<-h[[1]]
        t0<-unlist(lapply(seq(min(h),max(h)),
                          function(i)length(which(h==i))))
        p1<-heightHist(t0,which.max(t0)+min(h)-1)
        p2<-locHist(h,paste(m1,m2,sep="-"))
        print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
        print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=2))
    }
    else{
        if(is.list(n)) sequ<-n
        else sequ<-seq(n)
        for( k in sequ){
            t0<-unlist(lapply(seq(min(h[[k]]),max(h[[k]])),
                              function(i)length(which(h[[k]]==i))))
            p1<-heightHist(t0,which.max(t0)+min(h[[k]])-1)+ylab(name[k])
                tm2<-m2[k]
                tm1<-m1[k]
            p2<-locHist(h[[k]],paste(tm1,tm2,sep="-"))+ylab("")
            print(p1,vp=viewport(layout.pos.row=k,layout.pos.col=1))
            print(p2,vp=viewport(layout.pos.row=k,layout.pos.col=2))
        }
    }
}

#' @title Height Histrogram
#'  A plotting tool for histograms, Requires GGPLOT2
#' @export
heightHist<-function(t0,xlab="Histogram"){
    p<-ggplot(as.data.frame(t0),aes(x=t0))+geom_histogram(binwidth=1,xlab=xlab)
    if((max(t0)-min(t0))<20)
        p<-p+scale_x_continuous(breaks=seq(min(t0), max(t0),1))
    p<-p+stat_function(fun=function(x) length(t0)*2^(-x/sqrt(var(t0))),colour="blue" )
    p<-p+xlab(xlab)+ylim(c(0,length(which(t0==0))))
    p}

#' @title Location Histogram
#' A plotting tool for histograms, Requires GGPLOT2
#' @export
locHist<-function(t0,xlab="Histogram",limits=c(-32,32)){
    x<-seq(min(t0),max(t0))
    t1<-unlist(lapply(x,function(x)length(which(x==t0))))
    locHist2(t1,x,xlab,limits)
    #p<-qplot(x,t1,geom="step",ylab="frequency",xlab=xlab)+stat_function(fun=function(x){0})#+xlim(limits)
    #n<-4
    #m<-ceiling(log2(abs(limits[2]-limits[1])/2))
    #s<-seq(max(m-3,1),m)
    #b<-c(-sapply(s,function(x)2^x),sapply(s,function(x)2^x))
    #p+scale_x_continuous(breaks=b,limits=limits)
}

#' @export
locHist2<-function(t1,x,xlab="Histogram",limits=c(-32,32)){
    p<-qplot(x,t1,geom="step",ylab="frequency",xlab=xlab)+stat_function(fun=function(x){0})#+xlim(limits)
    n<-4
    m<-ceiling(log2(abs(limits[2]-limits[1])/2))
    s<-seq(max(m-3,1),m)
    b<-c(-sapply(s,function(x)2^x),sapply(s,function(x)2^x))
    p+scale_x_continuous(breaks=b,limits=limits)
}
