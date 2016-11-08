#' map zip list
#'
#' Takes a list of lists and groups the nth member of each list
#' into a new one.
#' @param ... an arbitraraly long set of equal length lists
#' @examples
#' mapziplist(list(1,2,3),list("A","B","C"))
#' # [[1]]
#' # [1] "1" "A"
#' # 
#' # [[2]]
#' # [1] "2" "B"
#' # 
#' # [[3]]
#' # [1] "3" "C"
._mapziplist<-function(...){
    inlist<-list(...)
    l<-length(inlist[[1]])
    zip<-function(i)
        unlist(lapply(inlist,"[[",i))
    lapply(seq(l),zip)
}


#' hg19Sort
#'
#' Reorders the chro factor in data to that which is outputed by BWA.
#' This allows for the sorting step to be skipped, however if another chromasome is used then a new order must be defined.
#' 
#' @param data [chr,start,end]
#' @return The same members of data but sorted on chr and start
#' 
._hg19Sort<-function(data){    
    neworder<-Filter(function(x) x %in% levels(data$chro), strsplit("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chrX chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr20 chrY chr19 chr22 chr21 chrM" ," ")[[1]])
    data$chro<-factor(data$chro,neworder)
    data<-data[with(data,order(chro,start)),]
    data
}


#' Quantile Normalization
#'
#' Applies quantile normalization to data in matrix form.
#' 
#' \itemize{
#' \item{1. Order each of the columns}
#' \item{2. Sum each orderd row and normalize by the width}
#' \item{3. Put the data back in order}
#' }
#' @param data matrix of data to be normalized
#' @return matrix of the same dimensions as data
#' @examples
#' matr<-cbind(rbind(5,2,3,4),rbind(4,1,4,2),rbind(3,4,6,8))
#' qn(matr)
._qn <-function(data){
    shape<-dim(data)
    sequence<-apply(data,2,order)
    reverseSequence<-unlist(apply(sequence,2,order))
    ranks<-apply(matrix(unlist(lapply(seq(shape[2]),function(i,x,y) x[y[,i],i],data,sequence)),ncol=shape[2]),1,sum)/shape[2]
    apply(reverseSequence,2,function(x) ranks[x])}


#' Normalize
#'
#' privides unit normalization for numeric x
#' @param x numeric input
#' @return numeric of length x
#' @examples
#' x<-rnorm(100,10,2)
#' y<-normalize(x)
#' print(c(mean(y),sd(y)))
._normalize<-function(x){
    norm<-function(x){
                      (x-mean(x))/(sqrt(var(x)))}
    if(is.null(dim(x)))
        norm(x)
    else
        apply(x,2,norm)
}


#' Merge Fun
#'
#' 
._mergeFun<-function(ma,swapFun=swapFun){
    newCols<-swapFun(colnames(ma))
    #print(newCols)
    unc<-unique(newCols)
    outL<-c()
    for(i in unc){
        pos<-which(i==newCols)
        #print(pos)
        #print(data.frame(c=colnames(ma)[pos],v=do.call(rbind,lapply(pos,function(x) sum(ma[,x])))))
        if(length(pos)>1)
            temp<-CCCA:::._orM(ma[,pos])
        else{
            temp<-ma[,pos]
            #print(sum(temp))
        }
        outL<-cbind(outL,temp)}
    
    colnames(outL)<-unc
    outL
}

#' orM
#'
#' applies the or (|) function accross a matrix recursivly
#' @seealso ~\code{\link{ors} \link{"|"}}
#' @param mat input matrix of logicals
#' @return a vector of logicals
#' @examples
#' a<-as.matrix(cbind(c(TRUE,TRUE,TRUE),
#'                     c(FALSE,TRUE,FALSE)))
#' result<-orM(a)
#' result==c(TRUE,TRUE,TRUE)
._orM<-function(mat){
    l<-dim(mat)[2]
    if(l==0)
        return(NULL)
    if(l==1)
        return(mat[,1])
    if(l==2)
        return(mat[,1] | mat[,2])
    if(l>2){
        ltemp<-cbind((mat[,1] | mat[,2]),mat[,3:l])
        return ( orM(ltemp))
    }
}


#' PCA tp Matrix Transformation
#' 
#' @param x the normalized input data and eigenvectors 
#' @param n normalizing function applied to the eigevectors
#' @return the dot product of the normalized data and eigenvectors
._pca2Matr<-function(x,n=pass){
    if(is.null(x$normData) | is.null(x$eigenVectors))
        stop("x needs variables normData and eigenVectors")
                                        # need to add tests for dimensions
    t(x$normData)%*%apply(x$eigenVectors,2,n)
}


._plotPCMatAux<-function(df,pcs,categories,colours,sdf=NULL,text=NULL,legend=NULL,label=NULL,blank=NULL){
    if(is.null(sdf))
       postext<-df
   else
       postext<-shiftCols(df$x,df$y,categories,sdf)
    p<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+geom_point(size=10,shape=20)+theme_bw()+ylab(pcs[2])+xlab(pcs[1])
    if(! is.null(blank))
    p<-p+ theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
    if(! is.null(label))
        p<-p+ylab("")+xlab("")
    if(! is.null(text))
        p<-p+geom_text(x=postext$x,y=postext$y,show_guide=F,size=5)
    if(! is.null(colours))
        p<-p+scale_color_manual(values=colours)
    if(is.null(legend))
        p<-p+theme(legend.position="none")+scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)
    p
}
