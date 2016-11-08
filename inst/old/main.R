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

