#' makeCCCA
#'
#' Generate a ChIP Component object from peaks, reads, and a list of
#' categories
#' @param dataSets list of files containing the raw reads in bed format
#' @param peakLists list of files containing the peaks reads in MACS'
#' xls format
#' @param categories the moniker to be attached to each of the data set files
#' @param pvalue the pvalue that the AFS peaks will be filtered by
#' @return a list of class ccca containing afs,udm, and prc
#' @examples 
#' dataSets<-sapply(c("raw_sample1.bed","raw_sample2.bed","raw_sample3.bed"),
#'   function(file){
#'     system.file("extdata", file, package = "CCCA")
#'   })
#' peakLists<-sapply(c("sample1.xls","sample2.xls","sample3.xls"),
#'   function(file){
#'     system.file("extdata", file, package = "CCCA")
#'   })
#' categories<-c("s1","s2","s3")
#' ccca<-makeCCCA(dataSets,peakLists,categories,20)
#' ## Find dimensions that seperate data sets of interest
#' plot(ccca,1,2)
#' plot(ccca,1,3)
#' ccca<-local({
#'   nm<-CCCA::normalize(ccca$prc)
#'   ccca<-addReg(ccca,"s1",nm[,1]<(mean(nm[,1])-3*sd(nm[,1])))
#'   ccca<-addReg(ccca,"s2",nm[,1]>(mean(nm[,1])+3*sd(nm[,1])))
#'   ccca<-addReg(ccca,"s3",nm[,2]>(mean(nm[,2])+3*sd(nm[,2])))
#'   ccca<-addReg(ccca,"s1.me",ccca$reg[,"s1"] & !(ccca$reg[,"s2"]
#'                | ccca$reg[,"s3"]))
#'   ccca<-addReg(ccca,"s2.me",ccca$reg[,"s2"] & !(ccca$reg[,"s1"]
#'                | ccca$reg[,"s3"]))
#'   ccca<-addReg(ccca,"s3.me",ccca$reg[,"s3"] & !(ccca$reg[,"s2"]
#'                | ccca$reg[,"s1"]))
#'   ccca
#' })
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome=BSgenome.Hsapiens.UCSC.hg19
#' addFasta(ccca,genome)
#' }
#' @export
makeCCCA<-function(dataSets,peakLists,categories,pvalue){
    ##dataSets<-c("raw_sample1.bed","raw_sample2.bed","raw_sample3.bed")
    ##peakList<-c("sample1.xls","sample2.xls","sample3.xls")
    ##categories<-c("s1","s2","s3")
    afs<-makeAFS(peakLists,categories,pvalue=pvalue)
    udm<-makeUDM(afs,dataSets)
    prc<-makePRC(udm)
    ret<-list(afs=afs,udm=udm,prc=prc,fasta=NULL,reg=NULL,
              categories=categories)
    attr(ret,"class")<-"ccca"
    ret
}

#' Load CCCA
#'
#' Generate the ccca object from premade AFS and UDM files.
#' See \link[CCCA]{makeCCCA}
#' @param peaks the AFS tab file
#' @param heights The UDM tab file
#' @param categories an array of categories
#' @return A ccca object
#' @examples
#' afsfiles<-system.file("extdata", "sample.afs", package = "CCCA")
#' udmfiles<-system.file("extdata", "sample.udm", package = "CCCA")
#' ccca<-loadCCCA(afsfiles,udmfiles,c("s1","s2","s3"))
#' @export
loadCCCA<-function(peaks,heights,categories){
    afs<-readAFS(peaks)
    udm<-readUDM(heights)
    prc<-makePRC(udm)
    ret<-list(afs=afs,udm=udm,prc=prc,fasta=NULL,reg=NULL,
              categories=categories)
    attr(ret,"class")<-"ccca"
    ret

}

#' @method plot ccca
#' @export
plot.ccca<-function(x,PC1,PC2,
                    swapFun=function(x){x},
                    swapCat=function(x){x},...){
    matr<-._pca2Matr(x$prc)
    df<-data.frame(x=matr[,PC1],y=matr[,PC2],categories=x$categories,
                   Conditions=swapCat(x$categories))
    ._plotPCMatAux(df,c(PC1,PC2),x$categories,NULL,swapCat,...)    
}

## #' @method print ccca
## #' @export
## print.ccca<-function(){
## }

#' @method addReg ccca
#' @export
addReg.ccca<-function(x, tag,logic,...){
    if(is.null(x$reg)){
        x$reg<-cbind(logic)
        colnames(x$reg)=tag
    }
    else{
        if(length(logic)!=dim(x$reg)[1])
            stop("Reg length not equal to region being added")
        logicp<-cbind(logic)
        colnames(logicp)<-tag
        x$reg<-cbind(x$reg,logicp)
    }
    x
}


#' @method contribution ccca
#' @export
contribution.ccca<-function(ccca,i,swapFun=function(string)string,
                            swapColour=NULL,...){
    PC<-ccca$prc$eigenVectors[,i]
    over<-as.data.frame(ccca$afs)
    ._stackedContrib(PC, "contrib2",
                     ._mergeFun(over[4:dim(over)[2]],swapFun),
                     swapFun=swapFun,colourOveride =swapColour,...)
}

#' @method addFasta ccca
#' @export
addFasta.ccca<-function(ccca,genome,width=200,...){
    requireNamespace("Biostrings")
    if (is.null(ccca$afs$chr) | is.null(ccca$afs$start)) 
        stop("addFasta env list must contain afs$chr afs$start and afs$end")
    if (!requireNamespace("Biostrings")) 
        stop(paste0("Must install the Biostrings package from",
                    " Bioconductor.\n",
                    "source(\"https://bioconductor.org/biocLite.R\");",
                    " biocLite(\"Biostrings\")"))

    start<- (ccca$afs$start + ccca$afs$end)/2-floor(width/2)
    ccca$fasta <- getSeq(genome,
                         ccca$afs$chr,
                         start = start,
                         width = width)
    ccca
}


