#' Shift from zero
#' 
#' Checks if the result of subtracting width from value is zero
#'
#' @param value integer or double which cannot be less than 0
#' @param width integer or double
#' @return a range value+width value-width where value-width >0
._shiftFromZero<-function(summit){
   testBedSE<-cbind(summit-350,summit+350)
   x<-which(testBedSE[,1]<0)
   testBedSE[x,2]=testBedSE[x,2]-testBedSE[x,1]
   testBedSE[x,1]=0
   testBedSE
}


#' sort data frame
#'
#' A small utility to functionalize the sorting of data frames
#' @export
#' @examples
#' dd<-data.frame(initial1=LETTERS[runif(100,1,24)],
#'                initial2=LETTERS[runif(100,1,24)],
#'                age=floor(runif(100,21,35)))
#' sortDataFrame(dd,"initial2","age")
._sortDataFrame<-function(dd,...)    
    dd[do.call(order,lapply(list(...),function(x) dd[x])),]


#' unify bed file
#'
#' unify bed file
._unifyBedFile<-function(dataset,overlapWidth, chr="chr",summit="summit",name="name"){
    chrcats<-levels(dataset[chr,])
    l<-nrow(dataset)
    t1<-dataset[1:l-1,]
    t2<-dataset[2:l,]
    b<-(t1[,chr]==t2[,chr] & abs(as.numeric(t1[,summit])-as.numeric(t2[,summit]))<overlapWidth)
    prePeaks<-unlist(lapply(which(b==FALSE),function(x){c(x,x+1)}))
    peaks<-t(matrix(c(1,prePeaks[1:length(prePeaks)-1]),nrow=2))
    data<-CCCA:::._unityOutput(peaks,
                as.integer(dataset[,chr]),
                as.integer(dataset[,summit]),
                as.integer(dataset[,name]))
    names(data)<-list(chr,summit,"matrix")
    data[[chr]]<-as.factor(data[[chr]])
    levels(data[[chr]])<-levels(dataset[,chr])
    colnames(data$matrix)<-levels(dataset[,name])
    od<-data.frame(chr=data[[chr]],summit=data[[summit]])
    od<-cbind(od,data$matrix)
    od
}


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

#' Plot PCA Matrix Auxilary
#'
#' 
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


#' genSwapFun
#' 
#' takes two equal length vectors and returns a function that maps
#' each input to its corresponding output. This can be used to alias
#' vectors if the mapping is known beforehand.
#' Note: the input vector must consist of strings. 
#' @param ouput The output mapping
#' @param input The input mapping
#' @return A function that maps each input to a single output
#' and returns a string
#' @examples
#' swapFun<-genSwapFun(c("a","b"),c(1,2))
#' swapFun("a")==1
#' swapFun("b")==2
#' swapFun("c")==NA
._genSwapFun<-function(input,output){
    # old name getSwapCats
    if(!length(output)==length(input))
        warning("output and input lengths are not equal")
    if(! is.character(input))
        warning("inputs must be a character, i.e. a vector of strings")
    names(output)<-input
    function(x){as.character(unlist(lapply(x,function(x) output[x])))}
}


#' Pass
#'
#' This function is used as a placeholder when a filter is required.
#' @param x variable
#' @return x unchanged
#' @examples
#' a<-1
#' pass(a) #> 1
._pass<-function(x) x

#' Create Zip
#'
#' takes a vector and retuns a list of pairs
#' @param x even length input vector
#' @return a list of pairs
#' @examples
#' createZip(c(1,2,3,4,5,6)) == list(c(1,2),c(3,4),c(5,6))
._createZip<-function(x){
    is.even<-function(x)
        (length(x) %% 2) ==0
    if(! is.even(x)){
        warning("length(x) is not even.")
        # x<-rbind(x,x) # This option created too many issues.
                        # Just let it fail.
    }
    Map(function(i,j)cbind(x[i],x[j]),seq(1,length(x),2),seq(2,length(x),2))
}

#' Split Zip
#'
#' Extracts the first and second values from a list of pairs.
#' @examples
#' splitZip(createZip(c(1,2,3,4,5,6))) == list(c(1,3,5),c(2,4,6))
._splitZip<-function(inList)
    list(sapply(inList,"[",1),sapply(inList,"[",2))

