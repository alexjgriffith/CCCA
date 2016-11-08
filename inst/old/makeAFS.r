#' Make AFS
#' 
#' Take a list of peaks and create a unified set.
#'
#' @param peakList
#' @param categories
#' @param format
#' @param pValue
#' @export
makeAFS<-function(peakList,categories,format="xls",pvalue=NULL){
    if(format!="xls" & !is.null(pvalue)){
        stop("to have a pvalue cut off the file format must be in MACS XLS output")
    }
    loadBedFun<-switch(format,
                       xls=readPeaksXLS,
                       bed=readPeaksBed)
    testBed<-unifyBedFile(
        sortDataFrame(do.call(rbind,Filter(function(x) ! is.null(x), lapply( mapziplist(peakList,categories),function(x) do.call(loadBedFun, as.list(x))))),"chr","summit"),700)
    #return(testBed)
    testBedSE<-shiftFromZero(testBed$summit)
    width<-dim(testBed)[2]
    retData<-data.frame(chr=testBed[,1], start=testBedSE[,1],end=testBedSE[,2])
    #retData<-data.frame(chr=testBed[,1], summit=testBed$summit)
    cbind(retData,testBed[,3:width])
