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
# qn
# pca
# peakPartitions
# applyPeakPartitions
# plotPCs

#' Quantile Normalization
#'
#' Applies quantile normalization to data in matrix form.
#' 
#' \itemize{
#' \item{1. Order each of the columns}
#' \item{2. Sum each orderd row and normalize by the width}
#' \item{3. Put the data back in order}
#' }
#' @param data matrix of data to be normalized
#' @return matrix of the same dimensions as data
#' @export
qn <-function(data){
    shape<-dim(data)
    sequence<-apply(data,2,order)
    reverseSequence<-unlist(apply(sequence,2,order))
    ranks<-apply(matrix(unlist(lapply(seq(shape[2]),function(i,x,y) x[y[,i],i],data,sequence)),ncol=shape[2]),1,sum)/shape[2]
    apply(reverseSequence,2,function(x) ranks[x])}


#' PCA
#'
#' Normalizes a data matrix and then preformes quantile normalization
#' using \code{prcomp} on the transform of the normalized data. The
#' eigenvectors are returned in a list with the normalized data
#' @param data The matrix to be analyzed
#' @param norm The normalization method
#' norm may either be a user provided function or one of the
#' following
#' \describe{
#' \item{rowSumOne}{ all rows sum to 1}
#' \item{colSumOne}{ all columns sum to 1}
#' \item{normRow}{all rows have mean 0 and sd 1}
#' \item{normCol}{all cols have mean 0 and sd 1}
#' \item{rowsVarOne}{ all rows have sd 1}
#' \item{colsVarOne}{ all cols have sd 1}
#' \item{qn}{ quantile normalization}
#' \item{none}{ no normalization}
#' }
#' @return list($normData,$eigenVectors)
#' @export
pca<-function(data,norm="qn"){
    print(head(data))
    if( is.function(norm))
       normData<-norm(data)
    else
        normData<-switch(norm,
               rowSumOne=t(apply(data,1, function(x) {x/sum(x)})),
               colSumOne=apply(data,2, function(x) x/sum(x)),
               normRow=t(apply(data,1, function(x) (x-mean(x))/var(x))),
               normCol=apply(data,2, function(x) (x-mean(x))/var(x)),
               rowsVarOne=t(apply(data,1, function(x) x/var(x))),
               colsVarOne=apply(data,2, function(x) x/var(x)),
               qn=qn(data),
               none=data,
               data)
    print(head(normData))
    prc<-prcomp(t(normData))$rotation
    list(normData=normData,eigenVectors=prc)
}

#' peak partitions
#'
#' Partitions a vector of weights based on its mean and standard deviation.
#' @param pcs the vectors to be assesed
#' @param l the number of vectors being assesed
#' @param tests the tests being applied to each vector
#' @param ns the number of stadard divations from the mean is considered
#' significant
#' @export
peakPartitions<-function(pcs,l=mWidth(pcs),tests=rep("gt",l),ns=rep(1,l)){
    tget<-function(x,i){
        if(l==1){
            x
        }else{
            x[,i]}
    }
    gt<-function(i){
        y<-tget(pcs,i)
        (y>(ns[i]*sqrt(var(y))+mean(y)))
    }
    lt<-function(i){
        y<-tget(pcs,i)
        (y<(-ns[i]*sqrt(var(y))+mean(y)))}
    test<-function(x)
        do.call(tests[x],as.list((x)))
    dta<-lapply(seq(l),test)
    if(l==1)
        return(unlist(dta))
    else
        do.call(and,list(dta))
}

#' apply peak partitions
#'
#' a wrapper for peak partitions that takes an inObj that contains all
#' of the information nessicary for the analysis
#'
#' @param prc the vectors
#' @param inObj list(list(which vector, function, ns) ...)
#' @return the regions seperated
#' @export
applyPeakPartitions<-function(prc,inObj){
    ret<-c()
    for (i in inObj)
        ret<-cbind(ret,peakPartitions(prc[,i[[1]]],test=i[[2]],ns=i[[3]]))
    return(ret)
}

#' plot PCs
#'
#' A tool used to plot the cross product of the normalized data
#' and the eigenvectors, it relies only on the standard library
#' ploting function
#' @param pcs the eiginvectors
#' @param pos a list indicating which two dimensions should be ploted
#' @param data the normalized data
#' @param cats the catagories, used to label the plot
#' @param lab the labeling information
#' @examples
#' heightFile<-"~/Dropbox/UTX-Alex/jan/combined_heights.bed"
#' catFile<-"~/Dropbox/UTX-Alex/jan/catagories"
#' cats<-as.character(unlist(read.table(catFile)))
#' score<-loadHeightFile(heightFile)$data
#' normData_eigen<-pca(score)
#' plotPCs(normData_eigen[[2]],c(1,3),normData_eigen[[1]],cats)
#' @export
plotPCs<-function(pcs,pos,data,cats,lab=c("xlabel","ylable","Title")){
    if("rotation" %in% names(pcs))
        x<-t(as.matrix(pcs$rotation)) %*% as.matrix(data)
    else
        x<-  t(as.matrix(pcs)) %*% as.matrix(data)
    d1<-data.frame(x[pos[1],])
    d2<-data.frame(x[pos[2],])
    plot(t(d1),t(d2),xlab=lab[1],ylab=lab[2])
    title(main=lab[3])
    text(x[pos[1],],x[pos[2],],labels=cats,cex=0.7,pos=3)}
