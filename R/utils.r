#!/usr/bin/env R
#
# This file is part of CCCA,
# http://github.com/alexjgriffith/CCCA/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com

#' ors
#'
#' Recursivly finds the union of a list of logicals
#' @export
ors<-function(inList){
    l<-length(inList)
    if(l==0)
        return(NULL)
    if(l==1)
        return(inList[[1]])
    if(l==2)
        return(inList[[1]] | inList[[2]])
    if(l>2){
        ltemp<-append(list(inList[[1]] | inList[[2]]),inList[3:l])
        return ( and(ltemp))
    }
}

#' Overlap
#'
#' Takes a matrix of logical inputs, each column representing a data set each
#' row representing a location on the genome,
#'
#' @param f n by m logical matrix
#' @return A m by m matrix of column overlaps
#' @export
#' @examples
#' lmatrix<-matrix(as.logical(c(1,1,1,1,1,1,0,0,1,0,1,0)),ncol=3)
#' overlap(lmatrix)
overlap<-function(f){
    permutations<-function(lin)
        as.vector(outer(lin,lin,function(x,y) paste(x,y,sep="")))
    overlapNorm<-function(x,f){
        a<-which(f[,x[1]]==1)
        b<-which(f[,x[2]]==1)
        length(intersect(a,b))/length(union(a,b))
    }
    width<-dim(f)[2]
    m<-matrix(apply(do.call(cbind,lapply(strsplit(permutations(seq(width)),split=""),as.numeric)),2, overlapNorm,f),ncol=width)
    colnames(m)<-colnames(f)
    rownames(m)<-colnames(f)
    m
}

#' map zip list
#'
#' Takes a list of lists and groups the nth member of each list
#' into a new one.
#' @param ... an arbitraraly long set of equal length lists
#' @export
#' @examples
#' mapziplist(list(1,2,3),list("A","B","C"))
#' # [[1]]
#' # [1] "1" "A"
#' # 
#' # [[2]]
#' # [1] "2" "B"
#' # 
#' # [[3]]
#' # [1] "3" "C"
mapziplist<-function(...){
    inlist<-list(...)
    l<-length(inlist[[1]])
    zip<-function(i)
        unlist(lapply(inlist,"[[",i))
    lapply(seq(l),zip)
}


#' sort data frame
#'
#' A small utility to functionalize the sorting of data frames
#' @export
#' @examples
#' dd<-data.frame(initial1=LETTERS[runif(100,1,24)],
#'                initial2=LETTERS[runif(100,1,24)],
#'                age=floor(runif(100,21,35)))
#' sortDataFrame(dd,"initial2","age")
sortDataFrame<-function(dd,...)    
    dd[do.call(order,lapply(list(...),function(x) dd[x])),]


# no longer used
or<-function(...){
    for(i in list(...))
        {
            if(! is.null(i))
                return(i)
        }
}

mWidth<-function(x)
    or(dim(x)[2],1)

mLength<-function(x)
    or(dim(x)[1],length(x))
    
and<-function(inList){
    l<-length(inList)
    if(l==0)
        return(NULL)
    if(l==1)
        return(inList[[1]])
    if(l==2)
        return(inList[[1]] & inList[[2]])
    if(l>2){
        ltemp<-append(list(inList[[1]] & inList[[2]]),inList[3:l])
        return ( and(ltemp))
    }
}
