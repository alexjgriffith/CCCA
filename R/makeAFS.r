#' Make AFS
#' 
#' Take a list of peaks and create a unified set.
#'
#' @param peakList A list of strings refering to file locations containing Peak data
#' @param categories The monikers these peak files will be refered to ass
#' @param format The format to be imported, either bed or xls
#' @param pvalue The pvalue cut off to be applied, only valid with xls format
#' @param width The min width summits must be from one another in order to cluster
#' @return The AFS generated from the input peak lists <chr><summit><catagory>
#' @examples
#' peakList<-sapply(c("sample1.xls","sample2.xls","sample3.xls"),
#'   function(file){
#'     system.file("extdata", file, package = "CCCA")
#'   })
#' categories<-c("s1","s2","s3")
#' afs<-makeAFS(peakList,categories,pvalue=20)
#' @export
makeAFS<-function(peakList,categories,format="xls",pvalue=NULL,width=700){
    if(format!="xls" & !is.null(pvalue)){
        stop("to have a pvalue cut off the file format must be in MACS XLS output")
    }
    if(! format %in% c("xls","bed") ){
        stop("format options are (xls|bed)")
    }
    loadBedFun<-switch(format,
                       xls=function(a,b=a){readPeaksXLS(a,b,pvalue)},
                       bed=readPeaksBed)
    ## Load each of th files
    files<-lapply( CCCA:::._mapziplist(peakList,categories),function(x) do.call(loadBedFun, as.list(x)))
    ## Unify the peaks based on width (default 700)
    testBed<-CCCA:::._unifyBedFile(
        CCCA:::._sortDataFrame(do.call(rbind,Filter(function(x) ! is.null(x), files)),"chr","summit"),width)
    ## Make sure none of the peaks have values less than 0
    testBedSE<-CCCA:::._shiftFromZero(testBed$summit)
    ## Return a data frame of form <chr><start><end><h1>...<hn>
    width<-dim(testBed)[2]
    retData<-data.frame(chr=testBed[,1], start=testBedSE[,1],end=testBedSE[,2])
    ret<-cbind(retData,testBed[,3:width])
    attr(ret,"class")<-"AFS"
    ret
}

## #' @method print AFS
## #' @export
## print.AFS<-function(x,...){
##     if(!require('data.table'))
##         print.default(x)
##     else
##         print(as.data.table(x))
## }

#' Read AFS
#' @export
readAFS<-function(fname){
    ##"~/Dropbox/Data/AFS/22_pvalue_20_control_combined.txt"
    ret<-CCCA:::._orderBed(read.table(fname,header=T))
    attr(ret,"class")<-c("AFS","data.frame")
    ret
}
