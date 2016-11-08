#' load Bed File
#' 
#' Calls file_length and read_bed from hg19Height.c and returns a 3 column bed file
#' @param file the location of the bed file of interest
#' @return data.frame $chro $start $end
#' @examples
#'  data<-loadBedFile(file="filelocation")  
._loadBedFile<-function(file){
    file<-normalizePath(file)
    if(!file.exists(file)){
        stop(paste0("Can't find file: ",file))
    }
    else{
        fileLength<-as.integer(.C("file_length",file,stringLength=integer(1))[[2]])
        results<-.C("read_bed",file,chro=character(fileLength),start=integer(fileLength),end=integer(fileLength))
        data.frame(chro=results$chro,start=results$start,end=results$end)
    }		
}	

#' Get Pile Up
#'
#' A minimal wrapper for the pileup function from hg19Height.c
#' returns a integer vector of computed read pileups
#' @param file The raw data file of interest
#' @param bed The preloaded bed infromation including bed$start and bed$end
#' @param chroms a list of chromosomes whos string values have been repalced with ranks
#' @param peakLength The length of the bed data provided
._getPileUp<-function(file,bed,chroms,peakLength){
    file<-normalizePath(file)
    if(!file.exists(file)){
        stop(paste0("Can't find file: ",file))
    }
    if(length(bed$start)!=length(chroms))
        stop("Chromosomes and start length differ")
    start<-as.integer(as.character(bed$start))
    end<-as.integer(as.character(bed$end))
    peaknum<-as.integer(peakLength)
    score<-as.integer(rep(0,peakLength))
    results<-.C("pileup",file,chrom=chroms,start=start,end=end,peaknum=peaknum,score=score)
    results$score
}
