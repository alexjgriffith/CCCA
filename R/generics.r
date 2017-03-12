#' Contribution Plot
#'
#' Generates the contribution histogram of the prc along a specified Principle Compoenent
#' @rdname contribution
#' @param prc The prc object to be ploted
#' @param i the PC to plot along
#' @param swapFun A function mapping the categoies of prc to another set of categories
#' @param swapColour A function that maps the results of swapFun to specifc colours
#' @return A ggplot2 object
#' @examples
#'
#' swapFun<-makeSwapfun("s1 sample1 s2 sample2 s3 sample3")
#' swapColour<-makeSwapfun("sample1 red sample2 blue sample3 orange")
#' contribution(CCCA_prc,swapFun,swapColour)
#' @export
contribution<-function(prc,i,swapFun=function(string)string,swapColour=NULL,...){
    UseMethod("contribution",prc)
}

## #' Contribution Plot
## clust<-function(x,...){
##     UseMethod("clust",x)
## }

#' Normalize Plot
#' @param x Vector or matrix to be normalized
#' @examples
#' r<-runif(100,0,100)
#' normalize(x)
#' b=NULL
#' b$eivenVector<-matrix(r,10)
#' class(b)<-"prc"
#' normalize(b)$eigenVector
#' @export
normalize<-function(x,...){
    UseMethod("normalize",x)
}

#' AddReg Plot
#'
#' @param x The object with a $reg value to be added
#' @param tag The moniker of the new region to be added
#' @param logic A logical list of equal length to other regions in x$reg
#' @export
addReg<-function(x, tag,logic,...){
    UseMethod("addReg",x)
}

#' AddReg Plot
#'
#' @param ccca The object used to generate the fasta strings
#' @param genome A BS.genome object 
#' @param width The width of fasta files, centered around the summit of each peak in ccca
#' @param ... Extension aguments
#' @return ccca with fasta values
#' @examples
#' 
#' @export
addFasta<-function(ccca,genome,width=200,...){
    UseMethod("addFasta",ccca)
}

#' @rdname contribution
#' @method contribution default
#' @export
contribution.default<-function(prc,i,...){
    stop("Currently only implemented for the prc class.")
}

## #' @rdname clust
## #' @method clust default
## #' @export
## clust.default<-function(x,...){
##     hclust(x,...)
## }

#' @rdname normalize
#' @method normalize default
#' @export
normalize.default<-function(x,...){
   ._normalize(x)
}

#' @rdname addReg
#' @method addReg default
#' @export
addReg.default<-function(x, tag,logic,...){
    stop("Currently only implemented for the ccca class.")
}

#' @rdname addFasta
#' @method addFasta default
#' @export
addFasta.default<-function(ccca,genome,width=200,...){
    stop("Currently only implemented for the ccca class.")
}
