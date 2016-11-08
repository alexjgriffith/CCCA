#' Make UDM
#'
#' Generates a pile up matrix from a unified set of peaks and a list of raw data sets
#' @param data A preloaded bed data.frame which includes slots $chro $start $end
#' @param rawdata a list of raw data files
#' @param n the number of nodes to use. If 0 then the parrallel package is not used
#' @param verboes if not null then print each file that is found
#' @param clust pass in a cluster defined outside of makeUDM
#' @return A matrix of class UDM containing the pile up counts under each peak in the AFS
#' @examples
#' library(paralell)
#' dataSets<-c("raw_sample1.bed","raw_sample2.bed","raw_sample3.bed")
#' cl<-makeForkcluster(3)
#' UDM<-makeUDM(afs,raw,n=1,clust=cl)
#' stopCluster(cl)
#' @export
makeUDM<-function(data,rawdata,n=0,verbose=NULL,clust=NULL){
    for(file in rawdata){
        if(! file.exists(file)){
            stop("Can't find file ",file,".")
            }
        if(!is.null(verbose))
            print(paste("# Raw data file ",file," was found.",sep=""))
    }
    peakLength<-length(data$chr)
    chroms<-as.character(data$chr)
    ## Parallell Implementation
    if(n>0){
        if(!require(parallel))
            stop("The parallel package must be loaded to run in paralell")
        ## if a cluster is passed in then make one now
        if(is.null(clust))
            cs<-makeForkCluster(n,renice=0)
        else
            cs<-clust
        ret<-matrix(unlist(parLapply(cs,rawdata,._getPileUp,data,chroms,peakLength)),nrow=peakLength)
        if(is.null(clust))
            stopCluster(cs)
    }
    ## Serial Implementation    
    else{
        ret<-matrix(unlist(lapply(rawdata,getPileUp,data,chroms,peakLength)),nrow=peakLength)
    }
    attr(ret,"class")<-"UDM"
    ret
}

#' @method print UDM
#' @export
print.UDM<-function(x,...){

}
