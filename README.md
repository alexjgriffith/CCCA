Cross Condition ChIP Analysis (CCCA)
==
Alexander Griffith, University of Ottawa
--

### Warning
**v0.1.0** is currently on **master**, this is the last stable production ready version of CCCA. **v0.1.0** does not pass R CMD check.

The current **develop** branch is a subset of the functionality of **v0.1.0** the notes presented here apply only to the develop branch, and not all of the functionality is fully tested. 

___

### Installing CCCA

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

### Building the UDM

```R
rawReadFiles<-list("cem.bam","jrk.bam","ery.bam","hsc.bam")
AFS<-read.table("AFS.tab",header=T)
## since the rawReadFiles can be very large the files are not completely
## loaded into memory in pileUp
UDM<-pileUp(AFS,rawReadFiles)

```

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

### Selecting Peaks

```R
load(UDM)
load(prc)

## Visual Inspection to Identify Dimensions

## Application of k-means clustering

## plot False Positives
prc

## Define rule for seperating

```

This file is part of CCCA, http://github.com/alexjgriffith/CCCA/, and is Copyright (C) University of Ottawa, 2015. It is Licensed under the three-clause BSD License; see LICENSE.txt.

