#' @export
genFasta<-function(fastaFile,fastaIndex,chro,start,end,width,saveFile){
    ff<-normalizePath(fastaFile)
    fi<-normalizePath(fastaIndex)
    sf<-normalizePath(saveFile)
    len<-length(start)
    fastadata<-saveFile
    temp<-.C("generateFasta",ff,fi,as.character(chro),as.integer(start),
       as.integer(end),as.integer(len),as.integer(width),
       fastadata=fastadata)
    #fastadata=temp$fastadata
    #summit=floor(start+end-width)/2
    #data.frame(chr=chro,start=summit,end=summit+width,fastadata)
    readDNAStringSet(saveFile)
 }



#' @export
genFastaIndex<-function(fastaFile,fastaIndex=NULL){
    fastaFile<-normalizePath(fastaFile)
    if(is.null(fastaIndex))
        fastaIndex<-tempfile()
    else
        fastaIndex<-normalizePath(fastaIndex)
    .C("generateIndex",fastaFile,fastaIndex)
    fastaIndex
}

#' @export
saveFastaIndex<-function(fastaIndex,newLoc)
    file.rename(fastaIndex,newLoc)

#' @export
bed2Fasta<-function(bed,fastaFile,width,fastaIndex=NULL,saveFile=NULL){
    if(is.null(saveFile)){
        saveFile<-tempfile()
        write.table("touch",saveFile);
    }
    ff<-normalizePath(fastaFile)
    if(is.null(fastaIndex))
        fi<-genFastaIndex(fastaFile,fastaIndex)
    else
        fi<-normalizePath(fastaIndex)
    genFasta(ff,fi,bed[,1],bed[,2],bed[,3],width,saveFile)    
}

#' @export
saveFasta<-function(test,file){
    locs<-paste(">",test$chr,":",test$start,"-",test$end,sep="")
    out<-rep("",length(locs)*2)
    out[seq(1,length(locs)*2-1,by=2)]<-locs
    out[seq(2,length(locs)*2,by=2)]<-as.character(unlist(test$fastadata))
    write.table(out,file,row.names=FALSE,
                col.names=FALSE,quote=FALSE)
}

