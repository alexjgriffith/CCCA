#' Read Peaks XLS format
#'
#' Read in peaks that MACs has output in its XLS format. Works with MACS 2.1 and 1.3.1
#' @param file The MACS xls output to be imported
#' @param name A string as a moniker for this data set
#' @param pvalue log10 pvalue used as a cut off, default 0
#' @return a data.frame <chr><summit><name> | NULL if file has no entries
#' @examples
#' # need to create sample.xls
#' filename<-system.file("extdata","sample1.xls", package = "CCCA")
#' readPeaksXLS(filename)
#' @export
readPeaksXLS<-function(file,name=file,pvalue=0){
    bedData<-read.table(file,header=TRUE,skip="#")
    bedData<-bedData[bedData$X.log10.pvalue>pvalue,]
    if(dim(bedData)[1]>0){
        ret<-data.frame(bedData$chr,bedData$abs_summit,name)
        colnames(ret)<-c("chr","summit","name")
    }
    else
        ret<-NULL    
    ret
}
