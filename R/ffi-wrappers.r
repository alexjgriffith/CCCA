#' Compare multiple ChIP data sets across experimental conditions
#'
#' @docType package
#' @name CCCA
#' @author "Alexander Griffith <griffitaj@@gmail.com>"
NULL

#' @useDynLib CCCA, .registration =TRUE, .fixes = c_
#' @import parallel
#' @import Biostrings
#' @import functional
#' @import graphics
#' @importFrom stats var
#' @importFrom utils read.table str write.table
#' @import ggplot2
#' @import abind
NULL


#' unity output
#'
#' Wrapper for the unityOutput C function.
#' 
#' @param peaks A 2xN matrix representing overlaping peak regions
#' @param intChr N length list of chomosomes represented as integers
#' @param intSummit The summit location 
#' @param intname The data set or name that the peak come froms in int form
._unityOutput<-function(peaks,intChr,intSummit,intname){
    nchr<-length(unique(intname))
    lpeaks<-dim(peaks)[1]
    retChr<-integer(lpeaks)
    retMatrix<-integer(lpeaks*nchr)
    retSummit<-integer(lpeaks)
    data<-.C(c_unityOutput,as.integer(intChr),
             as.integer(intSummit),
             as.integer(intname),
             as.integer(peaks[,1]),
             as.integer(peaks[,2]),
             as.integer(lpeaks),
             as.integer(nchr),chr=retChr,summit=retSummit,matrix=retMatrix)
    list(chro=data$chr,summit=data$summit,
         matrix=t(matrix(data$matrix,nrow=nchr)))
}


#' load Bed File
#' 
#' Calls file_length and read_bed from hg19Height.c and returns a 3 column
#' bed file
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
        fileLength<-file.info(file)$size
                    as.integer(.C(c_file_length,file,
                                  stringLength=integer(1))[[2]])

        results<-.C(c_read_bed,file,chro=character(fileLength),
                    start=integer(fileLength),end=integer(fileLength))
        data.frame(chro=results$chro,start=results$start,end=results$end)
    }
}

#' Get Pile Up
#'
#' A minimal wrapper for the pileup function from hg19Height.c
#' returns a integer vector of computed read pileups
#' @param file The raw data file of interest
#' @param bed The preloaded bed infromation including bed$start and bed$end
#' @param chroms a list of chromosomes whos string values have
#' been repalced with ranks
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
    results<-.C(c_pileup,file,chrom=chroms,start=start,end=end,
                peaknum=peaknum,score=score)
    results$score
}
