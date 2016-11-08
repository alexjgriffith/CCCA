#' Generate
#'
#' Used to build the AFS and UDM.
#' @param peakFiles
#' @param alignedData
#' @param macsCutOff
#' @param cl
#' @param categories
#' @param norm
#' @return a list containing <bed><over><heights><prc<normData><eigenVectors>>
#' @export
generate<-function(peakFiles,alignedData,macsCutOff=0,cl=NULL,
                   categories=peakFiles,norm=qn){
    if(!is.null(cl)){
        n<-cl
        cl<-makeForkCluster(cl)
    }
    else{
        n<-0
    }
    over<-makeAFS(peakFiles,categories,pValue=macsCutOff)
           #          ,c(c("chr","start","end"),categories))
    heights<-pileUp(over,alignedData,clust=cl,n=n)
    if(!is.null(cl))
        stopCluster(cl)
    prc<-pca(heights,norm=norm)
    list(over=over,bed=over[,1:3],heights=heights,prc=prc,
         categories=categories)    
}

#' Add Fasta
#'
#' Appends a Biostrings Object to ENV containing the base pair infromation
#' for each peak. The input env list is reqired to have a bed value that
#' contains a start and an chr. The mean of the start and end values is
#' extended by the width and used to find the base pair infromation of
#' interest
#' @param env A list containing the location of the peaks
#' @param genome the genome to which the locations apply
#' @param width the width of the peaks
#' @return env with an attached fasta parameter
#' @export
addFasta<-function(env,genome=BSgenome.Hsapiens.UCSC.hg19,width=300){
    if(is.null(env$bed$chr) | is.null(env$bed$start))
        stop("addFasta env list must contain bed$chr bed$start and bed$end")
    if(!require(Biostrings)  )
        stop("Must install the Biostrings package from Bioconductor.
source(\"https://bioconductor.org/biocLite.R\"); biocLite(\"Biostrings\")")
    env$fasta<-getSeq(genome,env$bed$chr,
                      start=(env$bed$start+env$bed$end)/2 -
                          floor(width/2),width=width)
    env
}

#' seperate axis
#'
#' Generates a function that can be used to isolate subsets of peaks using
#' the eigenVectors of the normalized data.
#' @param ... list of options list(name,fun,sd,pc)
#' @param inclusive if false the subsets are mutualy exclusive
#' @export
sepAxis<-function(...,inclusive=TRUE){
    inList<-list(...)
    sep<-function(i,env){
        with(i,{
            do.call(fun,list(normalize(env$prc$eigenVectors[,pc]),sd))
        })
    }
    fun<-function(env){
        do.call(cbind,addNames(lapply(inList,sep,env),
                               sapply(inList,"[[","name"),list=TRUE))
    }
    if(inclusive){
           return(fun)
       }
    else{        
       return( function(env){
            tot<-fun(env)
            n<-colnames(tot)
            
            d<-rep(FALSE,dim(tot)[2])
            addColnames(t(apply(tot,1,
                    function(x)
                        if(sum(x)==1){return(x)}
                        else{return(d)})),n)
        })           
    }
}

#' Get Motif Info
#'
#' Find the location of motifs for downstream analysis
#' @param env 
#' @param motifs The sequences being looked for
#' @return env with the locations of motifs locations of peaks and
#' set of locations that are unique
#' @export
getMotifInfo<-function(env,motifs){
    add<-with(env,{        
        eboxLoc<-sapply(motifs,grepMotifs,fasta)
        peakLoc<-lapply(seq(4,dim(over)[2]),Compose(sel2(over),as.logical,which))
        names(peakLoc)<-categories
        uniqueLocs<-apply(over[,4:dim(over)[2]],1,sum)==1
        regLoc<-lapply(seq(dim(reg)[2]),Compose(sel2(reg),which))
        return(list(eboxLoc=eboxLoc,peakLoc=peakLoc,uniqueLocs=uniqueLocs,
                    regLoc=regLoc))
    })
    append(env,add)
}
