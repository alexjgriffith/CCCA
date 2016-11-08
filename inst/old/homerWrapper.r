
                                        #!/usr/bin/env R
#
# This file is part of CCCA,
# http://github.com/alexjgriffith/CCCA/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com
#
# functions:
# homerWrapper
#

#' Homer Wrapper
#'
#' This is only a wrapper for the DENOVO functionality of homer2, homer must be installed on the system before this can be used.
#' 
#' @param sequences Biostring dnastring object
#' @param foreground The foreground subset of sequences, list
#' @param background The background subset of sequences, list
#' @param homerLocation The locatoin of the homer executable on your system
#' @param motifsFile The file to save the motifs in,
#' if none then the motifs are note saved
#' @param opts A string of additional parameters to be passed to homer
#' @return a pwm of the motifs found.
#' @export
homerWrapper<-function(sequences,foreground,background,homerLocation,motifsFile=FALSE,opts="-S 25 -len 6"){
    if (!is.character(motifsFile))
        motifsFile<-tempfile()
    treatmentFile<-tempfile()
    controlFile<-tempfile()    
    writeXStringSet(sequences[foreground],treatmentFile)
    writeXStringSet(sequences[background],controlFile)
    cmd<-paste(homerLocation, "denovo -i ",treatmentFile," -b ",controlFile,opts," > ",motifsFile,sep=" ")
    system(cmd)
    loadPWM(motifsFile,"homer")
}


homerWrapperKnown<-function(sequences,foreground,background,homerLocation,motifs,motifsFile=FALSE,opts=""){
    if (!is.character(motifsFile))
        motifsFile<-tempfile()
    treatmentFile<-tempfile()
    controlFile<-tempfile()
    motifInFile<-tempfile()
    
    ##write.table(motifs,file=motifInFile,quote=FALSE,col.names=FALSE,row.names=FALSE)
    cat(motifList2Homer(motifs),file=motifInFile)
    print(read.delim(motifInFile,sep="\n"))
    writeXStringSet(sequences[foreground],treatmentFile)
    writeXStringSet(sequences[background],controlFile)
    cmd<-paste(homerLocation, "known -i ",treatmentFile," -b ",controlFile," -m ",motifInFile," ",opts," > ",motifsFile,sep=" ")
    print(cmd)
    system(cmd)
    loadPWM(motifsFile,"homer")
}