#' Contribution Plot
#'
#' Generates the contribution histogram of the prc along a specified Principle Compoenent
#' @rdname contribution
#' @param ccca The ccca object to be ploted
#' @param i the PC to plot along
#' @param swapFun A function mapping the categoies of prc to another set of categories
#' @param swapColour A function that maps the results of swapFun to specifc colours
#' @param ... additional arguments
#' @return A ggplot2 object
#' @examples
#' swapFun<-makeSwapFun("s1 sample1 s2 sample2 s3 sample3")
#' swapColour<-makeSwapFun("sample1 red sample2 blue sample3 orange")
#' contribution(CCCA_prc,1,swapFun,swapColour,n=1,steps=10)
#' @export
contribution<-function(ccca,i,swapFun=function(string)string,swapColour=NULL,...){
    UseMethod("contribution",ccca)
}

## #' Contribution Plot
## clust<-function(x,...){
##     UseMethod("clust",x)
## }

#' Normalize Plot
#' @param x Vector or matrix to be normalized
#' @param ... additional arguments
#' @examples
#' r<-runif(100,0,100)
#' normalize(r)
#' b=NULL
#' b$eigenVector<-matrix(r,10)
#' class(b)<-"PRC"
#' normalize(b)
#' @export
normalize<-function(x,...){
    UseMethod("normalize",x)
}

#' AddReg Plot
#'
#' @param x The object with a $reg value to be added
#' @param tag The moniker of the new region to be added
#' @param logic A logical list of equal length to other regions in x$reg
#' @param ... additional arguments
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
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome=BSgenome.Hsapiens.UCSC.hg19
#' addFasta(CCCA_prc,genome)
#' }
#' @export
addFasta<-function(ccca,genome,width=200,...){
    UseMethod("addFasta",ccca)
}

#' @rdname contribution
#' @method contribution default
#' @export
contribution.default<-function(ccca,i,...){
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
