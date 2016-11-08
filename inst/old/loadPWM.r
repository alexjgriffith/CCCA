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
# loadPWM
#


#' Load PWM
#' 
#' This utiltiy is designed to load several PWM formats and create a uniform layout within R for analysis
#' @param fileLocation The location of the motif file
#' @param version the name of the program that outputed the file. Currently the only versions are homer and jaspar.
#' @return a unified pwm
#' @export
loadPWM<-function(fileLocation,version="homer")
{
    loadDataBuilder<-function(splitfun,header=FALSE,skip=1,id=">",
                              split=" ",region=1){
        function(fileLocation){
            data<-read.delim(fileLocation,header=header,skip=skip,sep="\n")
            out<-list()
            n<-0
            name<-c()
            info<-c()
            for(i in seq(length(t(data)))) {
                d<-strsplit(as.character(data[i,]),split)[[1]]
                l<-length(d)
                if(  id == strsplit(d[1],"")[[1]][region])
                    {
                        if (! n==0){out<-append(out,list(matrix(box,4)))}
                        box<-c()
                        name<-c(name,d[1])
                        info<-c(info,paste(d[2:l],sep="\t",collapse=""))
                        n<-n+1}
                else{
                    box<-c(box,splitfun(d,l))
                }}
            cbind(name=name,info=info,data=append(out,list(matrix(box,4))) )}}
    loadFunctions<-list(
        homer=c(
            skip=0,
            split="\t",
            splitfun=function(d,l){
                as.numeric(unlist(strsplit(as.character(d),"\t")))}),
        jaspar=c(skip=0,
            splitfun=function(d,l){
                na.omit(as.numeric(unlist(strsplit(unlist(d), "[^0-9]+"))))
        }))
    if (version %in% names(loadFunctions))
        do.call(loadDataBuilder,loadFunctions[[version]])(fileLocation)
    else
        warning(paste(version,"is not an option."))
}
