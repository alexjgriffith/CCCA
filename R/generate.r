makePeaks<-function(chrs,ranges,n){
    np<-floor(n/(length(chrs)))
    nf<-n-np*(length(chrs)-1)
    num<-c(nf,rep(np,(length(chrs)-1)))
    ret<-do.call(rbind,mapply(function(chr,max,n){
         data.frame(chr=chr,summit=ceiling(runif(n,25,max)))
    }, chrs,ranges,num,SIMPLIFY=FALSE))
    ret[order(ret$chr,ret$summit),]
}




makeReads<-function(peaks,chrs,ranges,s2n,n){
    rNoise<-floor((1-s2n)*n)
    rSignal<-ceiling(s2n*n)
    df1<-makePeaks(chrs,ranges,rNoise)
    rsp<-rSignal/length(peaks)
    df2<-do.call(rbind,mapply( function(chr,summit) {
        data.frame(chr=chr,summit=ceiling(runif(rsp, summit-200, summit+200)))
    },peaks$chr,peaks$summit ,SIMPLIFY=FALSE))
    df<-rbind(df1,df2)
    df$start=df$summit
    df$end=df$summit+50
    ret<-df[order(df$chr,df$start),]
    ret[c("chr","start","end")]
}


peaks<-makePeaks(c("chr1","chr2","chr3"),c(100000,100000,100000),100)
reads<-makeReads(peaks,c("chr1","chr2","chr3"),c(100000,100000,100000),0.5,100000)
