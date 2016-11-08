#' CCCA
#'
#' @param dataSets list of files containing the raw reads in bed format
#' @param peakLists list of files containing the peaks reads in MACS xls format
#' @param categories the moniker to be attached to each of the data set files
#' @return a list of class ccca containing afs,udm, and prc
#' @examples
#' dataSets<-c("raw_sample1.bed","raw_sample2.bed","raw_sample3.bed")
#' peakList<-c("sample1.xls","sample2.xls","sample3.xls")
#' categories<-c("s1","s2","s3")
#' ccca<-ccca(dataSets,peakLists,categories)
#' ## Find dimensions that seperate data sets of interest
#' plot(ccca,1,2)
#' plot(ccca,1,3)
#' ccca<-local({
#'   nm<-normalize(ccca$prc)
#'   ccca<-addReg(ccca,"s1",nm[,1]<(mean(nm[,1])-3*sd(nm[,1])))
#'   ccca<-addReg(ccca,"s2",nm[,1]>(mean(nm[,1])+3*sd(nm[,1])))
#'   ccca<-addReg(ccca,"s3",nm[,2]>(mean(nm[,2])+3*sd(nm[,2])))
#'   ccca<-addReg(ccca,"s1.me",ccca$reg[,"s1"] & !(ccca$reg[,"s2"] | ccca$reg[,"s3"]))
#'   ccca<-addReg(ccca,"s2.me",ccca$reg[,"s2"] & !(ccca$reg[,"s1"] | ccca$reg[,"s3"]))
#'   ccca<-addReg(ccca,"s3.me",ccca$reg[,"s3"] & !(ccca$reg[,"s2"] | ccca$reg[,"s1"]))
#'   ccca
#' })
#' ccca<-addFasta(ccca,genome)
#' s1Fasta<-ccca$fasta[ccca$reg[,"s1.me"]]
#' @export
ccca<-function(dataSets,peakLists,categories){
    dataSets<-c("raw_sample1.bed","raw_sample2.bed","raw_sample3.bed")
    peakList<-c("sample1.xls","sample2.xls","sample3.xls")
    categories<-c("s1","s2","s3")
    afs<-makeAFS(peakList,categories,pvalue=20)
    udm<-makeUDM(afs,dataSets)
    prc<-makePRC(udm)
    ret<-list(afs=afs,udm=udm,prc=prc,fasta=NULL,reg=NULL,categories=categories)
    attr(ret,"class")<-"ccca"
    ret
}

#' @method print ccca
#' @export
print.ccca<-function(){
}

#' @method addReg ccca
#' @export
addReg.ccca<-function(){
}


#' @method normalize ccca
#' @export
contributions.ccca<-function(prc,swapFun=function(string)string,swapColour=NULL,...){
    PC<-ccca$prc$eigenVectors
    over<-ccca$afs
    CCCA:::._stackedContrib(PC, "contrib2",CCCA:::._mergeFun(over[4:dim(over)[2]],swapFun),swapFun=swapFun,colourOveride =swapColour,...)
}

#' @method addFasta ccca
#' @export
addFasta.ccca<-function(ccca,genome,...){
    require('Biostrings')
}
