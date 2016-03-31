library(CCCA)
library(Biostrings)

#' Is Palandrome
#' Tests if biostring is a palandrome... takes IUPAC caracters
#' as input and returns a TRUE or FALSE
#' @param stringIn a string, limited to IUPAC characters
#' @return TRUE if palandrme FALSE if not
#' @examples
#' stringIn<-"CANNTG" # NN E-Box
#' ispalandrome(stringIn)
#' > TRUE
#' stringIn<-"GATTAG" 
#' ispalandrome(stringIn)
#' > FALSE
#' @export
ispalandrome<-function(stringIn){
    reverseString<-reverseIUPAC(stringIn)
    stringIn==reverseString
}

#' pariwise Combinations
#'
#' Takes a list and returns the pairwise combinations of all members
#' @param options list or cons
#' @return cbind(lista,listb)
#' @examples
#' characters<-c("A","C","G","T")
#' pairwiseCombs(characters)
#' @export
pairwiseCombs<-function(options){ 
    singleC<-do.call(rbind, as.list(rep(list(combn(options, 1)),2)))
    out1<-cbind(combn(options,2),singleC)
    out1[,order(out1[1,])] #orders the output
}

#' reverse IUPAC
#'
#' Takes one set of iupac characters and returns the reverse compliment
#' @param stringIn a string of iupac characters
#' @return returns a string equal in length to stringIn
#' @examples
#' motif<-"CANN"
#' reverseIUPAC("CANN")
#' # > NNTG
#' @export
reverseIUPAC<-function(stringIn)
    consenusIUPAC(compliment(IUPACtoBase(stringIn)))

#' permutations
#'
#' finds the permutations of all members in lin
#' @param lin (list input)
#' @return a vector of all permutations
#' @examples
#' permutations(c("a","c","g","t"))
#' @export
permutations<-function(lin)
    sort(as.vector(outer(lin,lin,function(x,y) paste(x,y,sep=""))))



genEboxPerms<-function(){
    options<-strsplit("ATGC","")[[1]]
    optcombs<-sort(apply(cbind(combn(options,2),rbind(combn(options,1),combn(options,1))),2,function(x) paste(x,collapse="")))
    unlist(lapply(permutations(options),function(x) paste("CA",x,"TG",sep="")))
}

genEboxCombs<-function(){
    options<-strsplit("ATGC","")[[1]]
    optcombs<-unique((function(x) apply(cbind(sapply(x,reverseIUPAC),x),1,function(x) sort(x)[1] ) )(permutations(options)))
    unlist(lapply(optcombs,function(x) paste("CA",x,"TG",sep="")))
}


grepMotifs<-function(stringIn,Sequences){
    motif<-IUPACtoBase(stringIn)
    compl<-compliment(motif)
    mac<-list(motif,compl)
    unlist(do.call(union, lapply(mac, grep,Sequences)))    
}

    
#' unique Tags
#' 
#' Takes a 4 x N data matrix representing the bed file
#' <chr><start><end><type1-type2-type3>
#' @param peaks
#' @param sep
#' @param column
#' @return a short cons containing the unique strings that make up the
#' peaks data frame level for the coulumn enterd
#' @examples
#' uniqueTags(peaks,column = 1) # to show the chromasomes
#' uniqueTags(peaks,column = 4) # to show the tags
#' @export
uniqueTags<-function(peaks,sep="-",column=4)
    unique(unlist(strsplit(levels(peaks[,column]),sep)))

#' uniqueOnGenome
#'
#' Identifies the peaks with only one tag
#' @export
uniqueOnGenome<-function(peaks,sep="-",column)
    sapply(strsplit(as.character(unlist(peaks[,column])),sep),
           function(x) length(x) ==1 )

motifFreqGen<-function(funIn){    
    function(motifs,types,lMotifs,lPeaks,...){
        fun<-Vectorize(funIn(...))
        b<-lMotifs
        a<-lPeaks
        freq<-outer(a,b,fun)
        colnames(freq)<-motifs
        rownames(freq)<-types
        t(freq)
        }
}

motifUniqueFrequency<-function(...){
    funIn<-function(uniqueLocs){function(x,y){
        uniqueMotLocSize<-length(intersect(intersect(x,y),which(uniqueLocs)))
        uniqueSize<-length(intersect(x,which(uniqueLocs) ))
        uniqueMotLocSize/uniqueSize
    }}
    motifFreqGen(funIn)(...)
}

motifAllFrequency<-function(...){
    funIn<-function(){function(x,y){
        length(intersect(x,y))/length(x )
    }}
    motifFreqGen(funIn)(...)
}



# FastaFile<-"~/Dropbox/UTX-Alex/jan/combined.fasta"
# peakFile<-"~/Dropbox/UTX-Alex/jan/combined_mock.bed"
# heightFile<-"~/Dropbox/UTX-Alex/jan/combined_heights.bed"
# catagoryFile<-"~/Dropbox/UTX-Alex/jan/catagories"
# peaks<-read.table(peakFile)
# Sequences <- readDNAStringSet(FastaFile, "fasta")
# temp<-read.table("~/Dropbox/UTX-Alex/jan/combined_heights.bed")
# data<-as.matrix(temp[,4:dim(temp)[2]])
# catagories<-as.character(unlist(read.table(catagoryFile)))
# colnames(data)<-catagories

#motifs<-genEboxCombs()
#types<-uniqueTags(peaks)    

#mloc<-sapply(motifs,grepMotifs,Sequences)
#ploc<-lapply(types,grep,peaks[,4])
#uniqueLocs<-uniqueOnGenome(peaks)



#allData<-data.frame(count=do.call(rbind,lapply(mloc,length)))
#motifUniqueFrequency(motifs,types,mloc,ploc,uniqueLocs)
#motifAllFrequency(motifs,types,mloc,ploc)
