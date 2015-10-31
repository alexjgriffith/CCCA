#!/usr/bin/env R
#
# This file is part of mulcal,
# http://github.com/alexjgriffith/CCCA/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com
#
# functions:
# IUPACtoBase
# compliment
# consenusIUPAC
# splitBlist
# motifString
# PWMtoCons
#

#' Turn IUPAC lables into basic neucleotides
#'
#' @param char A string of IUPAC characters
#' c("A","G","C","T","R","Y","S","W","K","M","B","D","H","V","N")
#' @param rl choice to return a value
#' \describe{
#' \item{\strong{TRUE}}{Each posible string is genearated}
#' \item{\strong{FALSE}}{Direct translation to brakets R=[AG]}}
#' @export
#' @examples
#' > IUPACtoBase("ARYS")
#' [1] "A[AG][CT][CG]"
IUPACtoBase<-function(char,rl=FALSE){
    IUPAC<-strsplit(char,"")[[1]]
    # Define standard IPAC Characters
    IUPACCharacters<-list("A","C","G","T",
                          c("A","G"),         #R
                          c("C","T"),         #Y
                          c("C","G"),         #S
                          c("A","T"),         #W
                          c("G","T"),         #K
                          c("A","C"),         #M
                          c("C","G","T"),     #B
                          c("A","G","T"),     #D
                          c("A","C","T"),     #H
                          c("A","C","G"),     #V
                          c("A","C","G","T")) #N
    names(IUPACCharacters)<-c("A","C","G","T","R","Y","S","W",
                              "K","M","B","D","H","V","N")
    # Test to see if all input charactersa are valid IUPAC
    if(any( ! IUPAC %in% names(IUPACCharacters))){
        stop("Input string contains non IUPAC characters.")
    }
    # Replace each character with a cons containing up to 4
    # nucleotides
    vals<-sapply(IUPAC,function(x){(IUPACCharacters[x])})
    # if rl is true base becomes a cons of all posible outcomes
    # from the iupac ex TW > c(TA,TT)
    if(rl){
        Base<-c("")
        for(i in vals){
            kit<-c()
            for(j in Base)
                kit<-c(kit,paste(j,i,sep=""))
            Base<-kit}}
    else{
        Base<-c("")
        for(i in vals)
            if(length(i)==1)
                Base<-paste(Base,i,sep="")
            else
                Base<-paste(Base,
                            "[",
                            do.call(paste, as.list(c(i,sep=""))),
                            "]",sep="")
    }
    return(Base)
}

#' reverse compliment
#'
#' Takes a string which can include base nucleaotides (AGCT) and
#' brakets ("[","]") and returns the reverse compliment including brakets 
#' @export
#' @examples
#' > compliment("[ACT]TA")
#' [1] "TA[AGT]"
#' > compliment(IUPACtoBase("STA"))
#' [1] "TA[CG]"
#' > consenusIUPAC(compliment(IUPACtoBase("STA")))
#' [1] "TAS"
compliment<-function(string){
    chars<-c("A","G","C","T","[","]")
    names(chars)<-c("T","C","G","A","]","[")
    paste(rev(sapply(strsplit(string,"")[[1]],
                     function(x){(chars[x])})),collapse="")
}

    
#' Consensus IUPAC
#'
#' The reverse of IUPACtoBase, takes base nucleotides and transforms
#' them them to their IUPAC equivelent
#' @export
#' @examples
#' > consenusIUPAC( IUPACtoBase("ARYS"))
#' [1] "ARYS"
consenusIUPAC<-function(mstring){
        IUPACCharacters<-list("A","C","G","T",c("A","G"),c("C","T"),c("C","G"),c("A","T"),c("G","T"),c("A","C"),c("C","G","T"),c("A","G","T"),c("A","C","T"),c("A","C","G"),c("A","C","G","T"))
        names(IUPACCharacters)<-c("A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N")
        IUPACc<-unlist(lapply(IUPACCharacters,paste,collapse=""))
        x<-splitBlist(mstring)
        paste(lapply(x,function(x)names(which(x==IUPACc))),collapse="")}

#' Split by List
#'
#' Applied to DNA nucleotides (ACGT). Returns a list with [] removed
#' @export
#' @examples
#' > consensus<-"AGCT[AGCT]G"
#' > splitBlist(consensus)
#'      buf buf buf buf buf    buf
#' [1,] "A" "G" "C" "T" "AGCT" "G"
splitBlist<-function(mstring){
    ret<-c()
    buf<-""
    data<-strsplit(mstring,"")[[1]]
    i<-1
    while (i <= length(data)){
        buf<-""
        if (data[i]=="["){
            while(data[i]!="]"){
                if (data[i] %in% c("A","C","G","T"))
                    buf<-paste(buf,data[i],sep="")
                i<-i+1}}
        else if (data[i] %in% c("A","C","G","T"))
            buf<-data[i]
        ret<-cbind(ret,buf)
        i<-i+1
    }
    ret
}

#' Motif String
#'
#' A scoring mechanism to convert PWM to strings
#' @export
#' @examples
#' > motifString(matrix(c(1,1,1,1,7,1,1,1),ncol=2))
#' [1] "[ACGT]A"
motifString<-function(x){
    paste(apply(x,2,function(x){
        DNA<-c("A","C","G","T")
        if(sum(x)==0)
            return (paste(c("[",sort(DNA),"]"),collapse=""))
        x<-x/sum(x)
        or<-order(x,decreasing = TRUE)
        t<-or
        if(x[or[1]]>=.6){return(DNA[or[1]])}
        else if(sum(x[or[1:2]])>=0.8){t<-or[1:2]}
        else if((sum(x[or[1:3]])>=0.95)){t<-or[1:3]}
        paste(c("[",sort(DNA[t]),"]"),collapse="")
    }
                ),collapse="")
}

#' PWM to Consensus Motif
#'
#' A caller for ConsesnusIUPAC, Recives PWM and returns a consusensus motif
#' The PWMs should be loaded through the loadPWM utility.
#' @export
PWMtoCons<-function(x)
    consenusIUPAC(motifString(x))


