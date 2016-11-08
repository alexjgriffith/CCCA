makePeaks<-function(chrs,ranges,n,macsRange=c(5,200)){
    np<-floor(n/(length(chrs)))
    nf<-n-np*(length(chrs)-1)
    num<-c(nf,rep(np,(length(chrs)-1)))
    ret<-do.call(rbind,mapply(function(chr,max,n){
         data.frame(chr=chr,abs_summit=ceiling(runif(n,25,max)),X.log10.pvalue=rnbinom(n,20,0.4))
    }, chrs,ranges,num,SIMPLIFY=FALSE))
    ret[order(ret$chr,ret$abs_summit),]
}




makeReads<-function(peaks,chrs,ranges,s2n,n){
    rNoise<-floor((1-s2n)*n)
    rSignal<-ceiling(s2n*n)
    dft<-makePeaks(chrs,ranges,rNoise)
    df1<-data.frame(chr=dft$chr,summit=dft$abs_summit)
    rsp<-ceiling(rSignal/length(peaks)*(peaks$X.log10.pvalue)/sum(peaks$X.log10.pvalue))
    print(rsp)
    print(length(rsp))
    print(dim(peaks)[1])
    reads<-mapply( function(chr,summit,n) {
        data.frame(chr=chr,summit=ceiling(rnorm(n, summit, 50)))
    },peaks$chr,peaks$abs_summit,rsp ,SIMPLIFY=FALSE)
    df2<-do.call(rbind,reads)
    df<-rbind(df1,df2)
    df$start=df$summit
    df$end=df$summit+50
    df$strand="-"
    ret<-df[order(df$chr,df$start),]
    ret[c("chr","start","end","strand")]
}

set.seed(400)
peaks1<-makePeaks(c("chr1","chr2","chr3"),c(100000,100000,100000),70)
shared1<-peaks1[sample(seq(dim(peaks1)[1]),30),]
peaks2<-makePeaks(c("chr1","chr2","chr3"),c(100000,100000,100000),70)
shared2<-peaks2[sample(seq(dim(peaks2)[1]),30),]
peaks3<-rbind(shared1,
              shared2,
              makePeaks(c("chr1","chr2","chr3"),c(100000,100000,100000),40))
peaks3<-peaks3[order(peaks3$chr,peaks3$abs_summit),]
peaks1p<-rbind(peaks1,shared2)
peaks1p<-peaks1p[order(peaks1p$chr,peaks1p$abs_summit),]
peaks2p<-rbind(peaks2,shared1)
peaks2p<-peaks2p[order(peaks2p$chr,peaks2p$abs_summit),]
reads1<-makeReads(peaks1p,c("chr1","chr2","chr3"),c(100000,100000,100000),0.3,100000)
reads2<-makeReads(peaks2p,c("chr1","chr2","chr3"),c(100000,100000,100000),0.8,100000)
reads3<-makeReads(peaks3,c("chr1","chr2","chr3"),c(100000,100000,100000),0.2,100000)
write.table(reads1,"raw_sample1.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
write.table(reads2,"raw_sample2.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
write.table(reads3,"raw_sample3.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
write.table(peaks1p,"sample1.xls",quote=FALSE,row.names=FALSE,sep="\t")
write.table(peaks2p,"sample2.xls",quote=FALSE,row.names=FALSE,sep="\t")
write.table(peaks3,"sample3.xls",quote=FALSE,row.names=FALSE,sep="\t")
