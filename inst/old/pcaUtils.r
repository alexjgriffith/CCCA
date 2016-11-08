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
# readCategories
# getSwapCats
# stringToSwap
# simpleSwapFun
# qstem



#' Read Categories
#'
#' read a single column file and return a vector of strings
#' @param filename the filename
#' @return a vector of strings 
#' @export
readCategories<-Compose(read.table,unlist,as.character)

#' genSwapFun
#' 
#' takes two equal length vectors and returns a function that maps
#' each input to its corresponding output. This can be used to alias
#' vectors if the mapping is known beforehand.
#' Note: the input vector must consist of strings. 
#' @param ouput The output mapping
#' @param input The input mapping
#' @return A function that maps each input to a single output
#' and returns a string
#' @examples
#' swapFun<-genSwapFun(c("a","b"),c(1,2))
#' swapFun("a")==1
#' swapFun("b")==2
#' swapFun("c")==NA
#' @export
genSwapFun<-function(input,output){
    # old name getSwapCats
    if(!length(output)==length(input))
        warning("output and input lengths are not equal")
    if(! is.character(input))
        warning("inputs must be a character, i.e. a vector of strings")
    names(output)<-input
    function(x){as.character(unlist(lapply(x,function(x) output[x])))}
}

#' string to swap fun
#'
#' takes a string and produces a swap fun from it. The string is
#' broken allong spaces.
#' @param x input string
#' @return a function mapping the even values of x to the odd
#' @examples
#' stringToSwap("a b c d")("b")=="a"
#' @export
stringToSwap<-function(x){
    do.call(genSwapFun,splitZip(createZip(strsplit(x," ")[[1]])))
}

#' simple swap fun
#'
#' this function is proforms the same role of stringToSwap, but the even and
#' odd strings are reversed.
#' @param char input string
#' @return a function mapping the odd values of x to the even
#' @examples
#' simpleSwapFun("a b c d")("a")=="b"
#' @export
simpleSwapFun<-function(char){
    stringToSwap(paste (rev(strsplit(char, " ")[[1]]),collapse=" "))
}

#' Order Bed
#'
#' sorts a data.frame based on the first and third columns
#' @param ret data frame
#' @return sorted data frame
#' @examples
#' orderBed(data.frame(c("a","a","c"),c(1,2,3),c(5,4,6)))
#' @export
orderBed<-function(ret)
    ret[order(as.character(ret[,1]),ret[,3]),]


#' Quick Stem Plot
#'
#' Abbreviates stem plotting for the analysis of composite motifs
#' @param y all locations before determining pileup
#' @param xlim This range that is to be ploted
#' @param ... remaining argument to stem and plot
#' @examples
#' x<-runif(1000)*100
#' qstem(x)
#' @export
qstem<-function(y,xlim=c(-32,32),...){
    yn<-y[y>xlim[1]&y<xlim[2]]
    yt<-getHeights(yn);    
    stem(do.call(seq,as.list(range(yn))),yt,xlim=xlim,...)
}
