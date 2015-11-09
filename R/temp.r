# this file is for temporary c-wrapping
# it is dependant on the inclusion of mulcal

#' @export
regionWrapper<-function(temp,chroms){
    width<-length(chroms)
    buffer<-as.ineger(rep(0,width))
    outMatrix<-integer((width+2)*(length-1))
    length<-length(temp[,1])
    .C("region",temp[,1],temp[,3],temp[,2],buffer,outMatrix,length,width)
}

#' @export
peakRegions<-function(data){
    chroms<-levels(data[,1])    
    subset<-data[,data[,1]==chroms[1]]
    temp<-rbind(cbind(subset[,2],subset[,"name"],0),
                cbind(subset[,3],subset[,"name"],1))
    regionWrapper(temp,chroms)
}


# utils
#' map zip list
#'  @export
mapziplist<-function(...){
    inlist<-list(...)
    l<-length(inlist[[1]])
    zip<-function(i)
        unlist(lapply(inlist,"[[",i))
    lapply(seq(l),zip)
}

# utils
#' sort data frame
#'  @export
sortDataFrame<-function(dd,...)    
    dd[do.call(order,lapply(list(...),function(x) dd[x])),]

# tags.r
#' unify bed file
#'  @export
unifyBedFile<-function(dataset,overlapWidth, chr="chr",summit="summit",name="name"){
    chrcats<-levels(dataset[chr,])
    l<-nrow(dataset)
    t1<-dataset[1:l-1,]
    t2<-dataset[2:l,]
    b<-(t1[,chr]==t2[,chr] & abs(as.numeric(t1[,summit])-as.numeric(t2[,summit]))<overlapWidth)
    prePeaks<-unlist(lapply(which(b==FALSE),function(x){c(x,x+1)}))
    peaks<-t(matrix(c(1,prePeaks[1:length(prePeaks)-1]),nrow=2))
    data<-unityOutput(peaks,
                as.integer(dataset[,chr]),
                as.integer(dataset[,summit]),
                as.integer(dataset[,name]))
    names(data)<-list(chr,summit,"matrix")
    data[[chr]]<-as.factor(data[[chr]])
    levels(data[[chr]])<-levels(dataset[,chr])
    colnames(data$matrix)<-levels(dataset[,name])
    od<-data.frame(chr=data[[chr]],summit=data[[summit]])
    od<-cbind(od,data$matrix)
    od
}

#' unity output
#'  @export
unityOutput<-function(peaks,intChr,intSummit,intname){
    nchr<-length(unique(intname))
    lpeaks<-dim(peaks)[1]
    retChr<-integer(lpeaks)
    retMatrix<-integer(lpeaks*nchr)
    retSummit<-integer(lpeaks)
    data<-.C("unityOutput",as.integer(intChr),as.integer(intSummit),as.integer(intname),as.integer(peaks[,1]),as.integer(peaks[,2]),as.integer(lpeaks),as.integer(nchr),chr=retChr,summit=retSummit,matrix=retMatrix)
    list(chro=data$chr,summit=data$summit,matrix=t(matrix(data$matrix,nrow=nchr)))
   }

#' testZero
#'
#' Checks if the result of subtracting width from value is zero
#'
#' @param value integer or double which cannot be less than 0
#' @param width integer or double
#' @return a range value+width value-width where value-width >0
#' @template authorTemplate
#' @export
testZero<-function(value,width){
    shift<-value-width
    if(shift>=0)
        return (c(value-width,value+width))
    else{
        return (c(value-width-shift,value+width-shift))}}

#' Shift from zero
#' @export
shiftFromZero<-function(summit){
   testBedSE<-cbind(summit-350,summit+350)
   x<-which(testBedSE[,1]<0)
   testBedSE[x,2]=testBedSE[x,2]-testBedSE[x,1]
   testBedSE[x,1]=0
   testBedSE
}


#' make afs
#' @export
makeAFS<-function(peakList,categories,format="xls",pValue=FALSE){
    if(format!="xls" & pValue){
        print("to have a pvalue cut off the file format must be in MACS XLS output")
        return(0)
    }
    readPeaksXLS<-function(file,name=file){
        bedData<-read.table(file,header=TRUE,skip="#")
        bedData<-bedData[bedData$X.log10.pvalue>pValue,]
        ret<-data.frame(bedData$chr,bedData$abs_summit,name)
        colnames(ret)<-c("chr","summit","name")
        ret
    }
    readPeaksBed<-function(file,name=file){
        # much faster, but far more error prone
        bedData<-loadBedFile(file)
        colnames(bedData)<-c("chr","start","end")
        summit<-bedData$start+floor((bedData$end-bedData$start)/2)
        cbind(bedData,name,summit)
    }
    loadBedFun<-switch(format,
                       xls=readPeaksXLS,
                       bed=readPeaksBed)
    testBed<-unifyBedFile(
        sortDataFrame(do.call(rbind,lapply( mapziplist(peakList,categories),function(x) do.call(loadBedFun, as.list(x)))),"chr","summit"),700)
    testBedSE<-shiftFromZero(testBed$summit)
    width<-dim(testBed)[2]
    retData<-data.frame(chr=testBed[,1], start=testBedSE[,1],end=testBedSE[,2])
    cbind(retData,testBed[,3:width])
    
    
}



#' testF
#' @export
testF<-function(){
# categories<-c("cd34_new","eryt","jurk","cem_1")
# categories<-as.character(unlist(read.table("~/Dropbox/UTX-Alex/jan/catagories")))
# peakXLSFiles<-paste("~/Dropbox/Data/august_peaks/",categories,"~combined_mock_peaks.xls",sep="")
#a<-makeAFS(peakXLSFiles,categories)
# rawReadFiles<-paste(categories,"/",categories,"_unique_nodupes.bed",sep="")
# print(dim(a))
#peakBedFiles<-paste("~/.peaktemp/",categories,"~combined",sep="")
#Rprof("inst/rprof/test-bed.r")
#a<-makeAFS(peakBedFiles,categories,"bed") 
#Rprof(NULL)

    library(CCCA)
    library(parallel)
setwd("/mnt/brand01-00/mbrand_analysis")
    
    for(pvalue in c(0,5,7.5,10,12.5,15)){
    categories<-c("cd34_new","eryt","jurk","cem_1")
    peakXLSFiles<-paste("peaks/october/",categories,"/combined_mock_peaks.xls",sep="")
    rawDataFiles<-paste("data_sets/",categories,"/",categories,"_unique_nodupes.bed",sep="")
    a<-makeAFS(peakXLSFiles,categories,pValue=pvalue)
    width<-dim(a)[2]
    colnames(a)<-c("chro",colnames(a)[2:width])
    a<-hg19Sort(a)
    a<-a[a$chro!="chrY"& a$chro!="chrM"& a$chro!="chrY",]
    score<-pileUp(a,rawDataFiles,4,TRUE)
    colnames(score)<-categories

    frontName<-function(x) paste("~/thesis-november/4x4-",pvalue,"-",x,".data",sep="")
    write.table(1-cor(score),frontName("cor"),quote=FALSE,sep="\t")
    write.table(1-cor(qn(score)),frontName("qn-cor"),quote=FALSE,sep="\t")
    write.table(1-overlap(a[,4:width]),frontName("ovr"),quote=FALSE,sep="\t")
    write.table(a,frontName("score-matrix"),quote=FALSE,row.names=FALSE,sep="\t")
    write.table(score,frontName("overlap-matrix"),quote=FALSE,row.names=FALSE,sep="\t")
    
    categories<-c("k562_1","eryt","jurk","cem_1")
    peakXLSFiles<-paste("peaks/october/",categories,"/combined_mock_peaks.xls",sep="")
    rawDataFiles<-paste("data_sets/",categories,"/",categories,"_unique_nodupes.bed",sep="")
    a<-makeAFS(peakXLSFiles,categories,pValue=pvalue)
    width<-dim(a)[2]
    colnames(a)<-c("chro",colnames(a)[2:width])
    a<-hg19Sort(a)
    a<-a[a$chro!="chrY"& a$chro!="chrM"& a$chro!="chrY",]
    score<-pileUp(a,rawDataFiles,4,TRUE)
    colnames(score)<-categories
    frontName<-function(x) paste("~/thesis-november/4x4-b-",pvalue,"-",x,".data",sep="")
    write.table(1-cor(score),frontName("cor"),quote=FALSE,sep="\t")
    write.table(1-cor(qn(score)),frontName("qn-cor"),quote=FALSE,sep="\t")
    write.table(1-overlap(a[,4:width]),frontName("ovr"),quote=FALSE,sep="\t")
    write.table(a,frontName("score-matrix"),quote=FALSE,row.names=FALSE,sep="\t")
    write.table(score,frontName("overlap-matrix"),quote=FALSE,row.names=FALSE,sep="\t")


    overlap<-abind(lapply(c(0,5,7.5,10,12.5,15),function(x)read.table(paste("4x4-",x,"-ovr.data",sep=""))),along=3)
    correlation<-abind(lapply(c(0,5,7.5,10,12.5,15),function(x)read.table(paste("4x4-",x,"-cor.data",sep=""))),along=3)

    plot(correlation[3,4,]-correlation[1,2,],type="l")
    
    plot(overlap[3,4,]-overlap[1,2,],type="l")
}
}

permutations<-function(lin)
    as.vector(outer(lin,lin,function(x,y) paste(x,y,sep="")))

overlapNorm<-function(x,f){
    a<-which(f[,x[1]]==1)
    b<-which(f[,x[2]]==1)
    length(intersect(a,b))/length(union(a,b))
}



overlap<-function(f){
    width<-dim(f)[2]
    m<-matrix(apply(do.call(cbind,lapply(strsplit(permutations(seq(width)),split=""),as.numeric)),2, overlapNorm,f),ncol=4)
    colnames(m)<-colnames(f)
    rownames(m)<-colnames(f)
    m
}

