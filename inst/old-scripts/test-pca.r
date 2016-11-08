library("Biostrings")
library("CCCA")

## Utils

cats<-as.character(unlist(read.table("~/Dropbox/UTX-Alex/jan/catagories")))
fasta<-readDNAStringSet("~/Dropbox/UTX-Alex/jan/combined.fasta")
heights<-loadHeightFile("~/Dropbox/UTX-Alex/jan/combined_heights.bed")


prc<-pca(heights$data,norm="qn")

inObj<-list(tall=list(1,"lt",1),
            eryt=list(1,"gt",1),
            ecfc=list(3,"gt",1),
            other=list(c(3,5),c("gt","gt"),c(1,1)),
            hspc=list(c(3,7),c("gt","lt"),c(1,1)),
            meka=list(c(3,7),c("gt","gt"),c(1,1)),
            diff=list(7,"gt",1))

regs<-applyPeakPartitions(prc,inObj)

