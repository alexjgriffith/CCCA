Cross Condition ChIP Analysis (CCCA)
==
Alexander Griffith, University of Ottawa
--

### Warning
**v0.1.0** is currently on **master**, this is the last stable production ready version of CCCA. **v0.1.0** does not pass R CMD check.

The current **develop** branch is a subset of the functionality of **v0.1.0** the notes presented here apply only to the develop branch, and not all of the functionality is fully tested. 

___
### To Update for v0.2.0

```R
## Functions that will be exported
## Ectract the functionality for these functions from the subset of
## CCCA v0.1.0
readPeaksXLS # copy
readPeaksBed # copy
makeAFS #copy
print.AFS #new
plot.AFS # new : plot pairwise correlation and pairwise overlap

makeUDM ## copy from pileUp
print.UDM # new :plot pairwise correlation 


## Class prc
## list normData eigenVecotrs cat2cont categories colours

prc # copy: pca (include categories for contributions)
plot.prc # new
print.prc # new
summary.prc # new
normalize.prc #new

contributions.prc # copy and simplify 
falseP.prc # new ish

clust.prc # new (this is going to be difficult)
summary.prcClust # new
print.prcClust # new

## Class ccca
## The class that hold the results of the analysis
ccca 

print.ccca # new
addReg.ccca # new to replace generate

addFasta.ccca # copy make sure to require('Biostrings') and class(genome) == BSGenome

```

---
### Installing CCCA

Install Dependencies

```R
library(parallel)

cs<-makeForkCluster(4)
packs<-list(
    "devtools",
    "functional",
	"abind",
	"ggplot2",
	"grid",
	"httr" )

parLapply(cs,packs,function(x) 
install.packages(x,type="source",repos='http://cran.us.r-project.org'))

   
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("BSgenome.Hsapiens.UCSC.hg19")


```

From within R
```R
library(devtools)
install_github("alexjgriffith/CCCA")
```

From the shell
```sh
## Clone the repository
clone https://github.com/alexjgriffith/CCCA.git .
cd CCCA

## Install 
R CMD INSTALL .
```

---
### Building the AFS

```R
library(CCCA)
## list the locations of the peak files
## xls implies output from MACS (v 1.* or 2.*)
peakFiles<-list("cem-peaks.xls","jrk-peaks.xls",
	"ery-peaks.xls","hsc-peaks.xls")
## The categories of the peak files	
categories<-c("cem","jrk","ery","hsc")
## Load the peak files
peaks<-lapply(peakFiles,readPeaksXLS)
AFS<-makeAFS(peakList,categories,pValue=20)

## Save the afs as a tab seprated file
write.table("AFS.tab",quote=FALSE,rowname=FALSE,sep"\t")

## Save the AFS as RData
save(AFS,file="AFS.RData")

```
---
### Building the UDM

```R
rawReadFiles<-list("cem.bam","jrk.bam","ery.bam","hsc.bam")
AFS<-read.table("AFS.tab",header=T)
## since the rawReadFiles can be very large the files are not completely
## loaded into memory in pileUp
UDM<-pileUp(AFS,rawReadFiles)

```
---
### Normalizing the UDM and Find PCs

```R
## Quantile normalize the PCS and apply SVD to find the PCA
load(UDM.RData)
colnames(UDM)
prc<-pca(UDM,norm=qn)

## values used in naming and colour coding plots
prc$cat2cont<-makeSwapFun("hsc HSPC jrk TALL cem TALL ery ERYT")
prc$colours<-makeSwapFun("HSPC orange TALL blue ERYT red")

## Print Summary of principle components
prc

## Plot the pairwise combinations of the first 4 principle components
plot(prc)

## Plot the pairwise combinations of the first 2 principle components
plot(prc,PC1=1,PC2=2)

## Plot the colour coded histograms 
contributions(prc,PC1=1,cat2cont,colours)

```
---
### Selecting Peaks

```R
load(UDM)
load(prc)
load(AFS)

## Application of k-means clustering

clust<-kclust(prc,3)
summary(clust)

## Visual Inspection to Identify Dimensions
plot(prc,sd=2,clusters=clust)

## How many peaks will be selected using a cut off of 2
summary(prc,sd=2)

## plot False Positives
falseP<-falsePos(prc,sd=2,clusters=clust)

plot(falseP)

## Define rule for seperating
ccca<-ccca(prc=prc,udm=UDM,afs=afs)

ccca<-addReg(ccca,"ERYT",pc=1,sd=-2,fun="<")
ccca<-addReg(ccca,"TALL",pc=1,sd=2,fun=">")
ccca<-addReg(ccca,"HSPC",pc=2,sd=-2,fun="<")

summary(ccca)

save(ccca, file="ccca.RData")
```

---
### Functional Analysis

For functional analysis there are four complimentary packages for CCCA
	1. HomerWrapper
	2. STAMPWrapper
	3. mT1
	4. RGREAT

I also make use of the DAVID web service API R bindings. 

#### First get the necisarry pacakges
```R
library(devtools)
## make sure the dependencies are installed

install_github("alexjgriffith/HomerWrapper")
install_github("alexjgriffith/STAMPWrapper")
install_github("alexjgriffith/mT1")
install_github("alexjgriffith/RGREAT")
package.install("RDAVIDWebService")
```

#### Motif Denovo Search 
```R
library(HomerWrapper)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
load(ccca.RData)

## For large files this can be memory intesnsive
ccca<-addFasta(ccca,BSgenome.Hsapiens.UCSC.hg19)

## replace this with the location of you homer2 install
homerLocation="/usr/bin/homer2"

homerWrapper(ccca$fasta, ccca$reg[,"ERYT"],ccca$reg[,"CENTER"] , homerLocation, motifsFile = "eryt_6.pwm", opts = "-S 25 -len 6");
homerWrapper(ccca$fasta, ccca$reg[,"HSPC"],ccca$reg[,"CENTER"] , homerLocation, motifsFile = "hspc_6.pwm", opts = "-S 25 -len 6");
homerWrapper(ccca$fasta, ccca$reg[,"TALL"],ccca$reg[,"CENTER"] , homerLocation, motifsFile = "tall_6.pwm", opts = "-S 25 -len 6");



```

#### Motif Annotation
```
library(HomerWrapper)
library(StampWrapper)
homer<-list(eryt=loadPWM(eryt_6.pwm),
	hspc=loadPWM(hspc_6.pwm),
	tall=loadPWM(tall_6.pwm))
	

stamp<-lapply(homer,function(h){
	motifs<-gsub(">","",h[,1])
	pwms<-h[,2]
	names(pwms)<-motifs
	stampTRANSFAC<-stampWrapper(pwms,"TRANSFAC_Fams")
	stampJASPAR<-stampWrapper(pwm,"JASPAR_Fams")
	list(transfac_6=stampTRANSFAC,jaspar_6=stampJASPAR)
})


db<-buildDb(homer,stamp)

write.table(db,file="stamp_db.txt",quote=FALSE,rownames=FALSE)

```

#### Motif Clustering
```R
library(HomerWrapper)
library(StampWrapper)

## this step is more important if you are looking for motifs of different 
## lengths

db<-read.table("stamp_db.txt",header=T)
cluster(motifs,


```

#### Composite Motifs
```R
library(mT1)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

fasta<-mT1_fasta
motifs<-unique(c("CANNTG","GATAA", mT1_jaspar[1:100]))

obj<-mT1(fasta,motifs)


```
#### Gene Assosiation
```R

library(RGREAT)
load("RefSeq.hg19.csv")

```


#### Gene Ontology 

---
This file is part of CCCA, http://github.com/alexjgriffith/CCCA/, and is Copyright (C) University of Ottawa, 2015. It is Licensed under the three-clause BSD License; see LICENSE.txt.

