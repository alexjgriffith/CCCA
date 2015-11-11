#!/usr/bin/env R
#
# This file is part of CCCA,
# http://github.com/alexjgriffith/CCCA/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com

#' loadBedFile
#' 
#' Calls file_length and read_bed from hg19Height.c and returns a 3 column bed file
#' @param file the location of the bed file of interest
#' @return data.frame $chro $start $end
#' @export
#' @template authorTemplate
#' @examples
#'  data<-loadBedFile(file="filelocation")  
loadBedFile<-function(file){
    file<-normalizePath(file)
    if(!file.exists(file)){
        sprintf("Can't find file %s.",file)}
    else{
        fileLength<-as.integer(.C("file_length",file,stringLength=integer(1))[[2]])
        results<-.C("read_bed",file,chro=character(fileLength),start=integer(fileLength),end=integer(fileLength))
        data.frame(chro=results$chro,start=results$start,end=results$end)
    }}

#' convertChroms
#' 
#' Calls valueChromosome from hg19Height.c which converts a chromosme rank to a string
#' @seealso For the reverse see \code{\link{rankChroms}}
#' @param lis a list of chromosmes to
#' @template authorTemplate
#' @export
convertChroms<-function(lis){
    n<-length(lis)
    .C("valueChromosomes",character(n),as.integer(n), lis)[[1]]}

#' @export
getChroms<-function(file,n=100){
    file<-as.character(normalizePath(file))
    if(!file.exists(file)){        
        stop(paste("Can't find file",file))
    }
    else{        
        return(Filter(function(x) x!="", .C("getChroms",file,choms=character(n))$choms))}
}

#' rankChroms
#' 
#' Calls rankChromosome from hg19Height.c which converts a chromosme string to its equivelent rank
#' @seealso For the reverse see \code{\link{convertChroms}}
#' @template authorTemplate
#' @export
rankChroms<-function(lis){
    n<-length(lis)
    .C("rankChromosomes",as.character(lis),as.integer(n), integer(n))[[3]]}

#' getPileUp
#'
#' A minimal wrapper for the pileup function from hg19Height.c
#' returns a integer vector of computed read pileups
#' @param file The raw data file of interest
#' @param bed The preloaded bed infromation including bed$start and bed$end
#' @param chroms a list of chromosomes whos string values have been repalced with ranks
#' @param peakLength The length of the bed data provided
#' @export
getPileUp<-function(file,bed,chroms,peakLength){
    start<-as.integer(as.character(bed$start))
    end<-as.integer(as.character(bed$end))
    peaknum<-as.integer(peakLength)
    score<-integer(peakLength)
    print(sapply(list(start,end,score,chroms),length))
    print(file)
    results<-.C("pileup",file,chrom=chroms,start=start,end=end,peaknum=peaknum,score=score)
    results$score
    #NULL
}

#' PileUp
#'
#' Generates a pile up matrix from a unified set of peaks and a list of raw data sets
#' @template authorTemplate
#' @param data A preloaded bed data.frame which includes slots $chro $start $end
#' @param rawdata a list of raw data files
#' @param n the number of nodes to use. If 0 then the parrallel package is not used
#' @examples
#' # Initialize catagories, files containing raw data 
#' cats<-read.table("/home/griffita/Dropbox/UTX-Alex/jan/catagories")
#' prefix<-"/mnt/brand01-00/mbrand_analysis/data_sets/"
#' suffix<-"_sorted.bed"
#' rawdata<-apply(cats,1,function(x){paste(prefix,x,"/",x,suffix,sep="")})
#' # Apply pileUp to peaks found using MACS
#' data<-hg19Sort(loadBedFile(file))
#' score<-pileUp(data,rawdata,n=22)
#' # cluster the data sets based on the read hights of the peaks
#' temp<-cor(score)
#' rownames(temp)<-t(cats)
#' colnames(temp)<-t(cats)
#' pdf("test.pdf")
#' plot(hclust(dist(temp)),hang=-1)
#' # Rather than using peaks we can do global analyis
#' # this relies on breaking the genome into bins
#' # human.hg19.genome provides the start and stop point for each genome
#' # bin size desides how large the regions should be
#' binSize=10000
#' regions<-do.call(rbind,
#' apply(read.table("/data/binaries/BEDTools/genomes/human.hg19.genome")[1:24,],1,
#' function(x,step) {y<-seq(1,as.numeric(x[2]),step);
#' cbind(as.character(x[1]),as.character(y),as.character(y+step))} ,binSize))
#'   data<-hg19Sort(data.frame(chro=regions[,1],
#'                             start=as.integer(regions[,2]),
#'                             end=as.integer(regions[,3])))
#'   score<-pileUp(data,rawdata,n=22)
#' @export
pileUp<-function(data,rawdata,n=0,verbose=FALSE){
    for(file in rawdata){
        if(!file.exists(file)){
            sprintf("Can't find file %s.",file)
            return}
        if(verbose)
            print(paste("# Raw data file ",file," was found.",sep=""))
    }
    peakLength<-length(data$chro)

    chroms<-rankChroms(data$chro)

    if(n>0){
        cs<-makeForkCluster(n,renice=0)
        ret<-matrix(unlist(parLapply(cs,rawdata,getPileUp,data,chroms,peakLength)),nrow=peakLength)
        stopCluster(cs)
    }else{
    ret<-matrix(unlist(lapply(rawdata,getPileUp,data,chroms,peakLength)),nrow=peakLength)}
    ret}

#' hg19Sort
#'
#' Reorders the chro factor in data to that which is outputed by BWA.
#' This allows for the sorting step to be skipped, however if another chromasome is used then a new order must be defined.
#' 
#' @param data [chr,start,end]
#' @return The same members of data but sorted on chr and start
#' @template authorTemplate
#' @export
hg19Sort<-function(data){    
    neworder<-Filter(function(x) x %in% levels(data$chro), strsplit("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chrX chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr20 chrY chr19 chr22 chr21 chrM" ," ")[[1]])
    data$chro<-factor(data$chro,neworder)
    data<-data[with(data,order(chro,start)),]
    data
}

#' Chrom Genome Sort
#'
#' Returns a function that reorders the peak data based on the chromosome order provided
#' @param chromOrder A list of chromosomes to be found in the data
#' @examples
#' sortedChroms<-strsplit("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chrX
#' chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17
#' chr18 chr20 chrY chr19 chr22 chr21 chrM" ," ")[[1]]
#' hg19Sort<-chromGenomeSort(sortedChroms)
#' bedData<-loadBedData("test.bed")
#' sortedBedData<-hg19Sort(bedData)
#' @export
chromGenomeSort<-function(chromOrder){
    rf<-function(data){
        ls<-levels(data$chro)
        neworder<-Filter(function(x) {x %in% ls}, chromOrder)
        data$chro<-factor(data$chro,neworder)
        data<-data[with(data,order(chro,start)),]
        data                       
        }
   return(rf)
}
    
#' getData
#'
#' loads data from file returns a sorted set with a height matrix
#' @export
#' @template authorTemplate
#' @param file peak file
#' @param rawdata a list of raw data files
#' @param n the number of nodes used default=0
getData<-function(file,rawdata,n=0){
    data<-hg19Sort(loadBedFile(file))
    score<-pileUp(data,rawdata,n=n)
    c(data=data,score=score)}
