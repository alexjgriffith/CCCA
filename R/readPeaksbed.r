#' Read Peaks BED format
#'
#' Read in peaks that MACs has output in its bed format. Works with MACS 2.1 and 1.3.1.
#' Cannot filter peaks based on macs score.
#' @param file The MACS xls output to be imported
#' @param name A string as a moniker for this data set
#' @return a data.frame <chr><summit><name> | NULL if file has no entries
#' @examples
#' # need to create sample.xls
#' filename<-system.file("extdata","raw_sample1.bed", package = "CCCA")
#' readPeaksBed(filename)
#' @export
readPeaksBed<-function(file,name=file){ 
    bedData<-._loadBedFile(file)
    colnames(bedData)<-c("chr","start","end")
    summit<-bedData$start+floor((bedData$end-bedData$start)/2)
    data.frame(chr=bedData$chr,summit=summit,name=name)
}
