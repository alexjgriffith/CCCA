setwd("~/Masters/CCCA")
install()

library(CCCA)

data<-loadBedFile("~/.peaktemp/tall_p1~combined")

data<-read.table("~/Dropbox/Data/22x22-pvalue=0.matrix",header=T)
bed<-data[,1:3]
categories<-colnames(data)[4:dim(data)[2]]
heights<-read.table("~/Dropbox/Data/22x22-pvalue=0_heights.matrix",header=T)
colnames(heights)<-categories
#test<-bed2Fasta(bed,"~/hg19.fasta",300,fastaIndex="~/hg19.index")

pcs<-pca(heights)

Erythroid<-peakPartitions(pcs$eigenVectors[,1],1,"gt",1)
TALL<-peakPartitions(pcs$eigenVectors[,1],1,"lt",1)


homerWrapper(test,Erythroid, !Erythroid,"~/Masters/mulcal/inst/lib/homer-4.7/bin/homer2", "inst/exdata/pv0-pc1-gt1-6.pwm")


motifs<-loadPWM("inst/exdata/eryt_vs_not_eryt.pwm")


library(functional)

sapply(motifs[,3],Compose(motifString,consenusIUPAC))


