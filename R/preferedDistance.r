#' find motif distance
#'
#' Identifies the indicies of the peaks that have two motifs at a prefered
#' distance from one another. Note that it only takes the smallest distance
#' into account when searching so this function is not suitible when the
#' expected distance between two motifs is great.
#' @param fasta The fasta sequences from Biostrings
#' @param motifa a string of IUPAC characters
#' @param motifb a string of IUPAC characters
#' @param loc a list of logicals equal in length to the fasta variable
#' @param ... the distances at which peaks will be returned
#' @return a list of indicies for peaks that have the motifs at prefered
#' distances
#' @export
findMotifDist<-function(fasta,motifa,motifb,loc,...){
    selectRegionsByDistance<-function(motif1,motif2,locations,fasta){
        m1<-gregexpr(motif1,fasta[locations],ignore.case=TRUE)
        m2<-gregexpr(motif2,fasta[locations],ignore.case=TRUE)
        c1<-gregexpr(compliment(motif1),fasta[locations],ignore.case=TRUE)
        c2<-gregexpr(compliment(motif2),fasta[locations],ignore.case=TRUE)
        findDmin<-function(a,b){
            c<-NA
            if(a[1]==-1 & b[1]!=-1){
                c<-NA#b[which.min(abs(b))]
            }
            else if (a[1]!=-1 & b[1]==-1){
                c<-NA#a[which.min(abs(a))]
            }
            else if (a[1]!=-1 & b[1]!=-1){
                c<-outer(a,b,"-")
               # NOTE:: by taking min we don't look at all possible distances
                c<-c[which.min(abs(c))]
            }
            else{
                c<-NA
            }
            c}
        cd<-(mapply(findDmin,c1,c2)+1)*-1
        md<-mapply(findDmin,m1,m2)
        a<-apply(cbind(md,cd),1,function(x)
            if(length(which(is.na(x)))>1){NA}
            else{min(x,na.rm=TRUE)})
        a
    }
    getMinDistance<-function(fasta,motifa,motifb,loc){
        a<-cbind(which(loc),
                 selectRegionsByDistance(IUPACtoBase(motifa),
                                         IUPACtoBase(motifb),
                                         loc,env$fasta))
        b<-a[!is.na(a[,2]),]
        b    
    }    
    selMotifReg<-function(minDistance,...){
        getMin<-function(loc,minDistance)
            minDistance[which(minDistance[,2]==loc),1]
        sort(unique(unlist(lapply(list(...),getMin,minDistance))))
    }    
    a<-getMinDistance(fasta,motifa,motifb,loc)
    selMotifReg(a,...)
}


