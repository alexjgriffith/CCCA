#!/usr/bin/env R
#
# This file is part of CCCA,
# http://github.com/alexjgriffith/CCCA/, 
# and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
# the three-clause BSD License; see LICENSE.txt.
# Author : Alexander Griffith
# Contact: griffitaj@gmail.com

# lzip
# pass
# createZip
# splitZip
# addRownames
# addColnames
# addNames
# normalize
# makeLogic
# stem
# ors
# orM
# unionN
# intersectN
# overlap
# mapziplist
# sortDataFrame
# mor
# mWidth
# mLength
# modulous


#' List Zip
#'
#' Zips 0-N lists
#' @examples
#' result<-lzip(list(1,2,3),list("a","b","c"))
#' result==list(list(1,"a"),list(2,"b"),list(3,"c"))
lzip<-function(...){
    apply(mapply(function(...)list(...),...),2,as.list)
}

#' Pass
#'
#' This function is used as a placeholder when a filter is required.
#' @param x variable
#' @return x unchanged
#' @examples
#' a<-1
#' pass(a) #> 1
#' @export
pass<-function(x) x

#' Create Zip
#'
#' takes a vector and retuns a list of pairs
#' @param x even length input vector
#' @return a list of pairs
#' @examples
#' createZip(c(1,2,3,4,5,6)) == list(c(1,2),c(3,4),c(5,6))
createZip<-function(x){
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
splitZip<-function(inList)
    list(sapply(inList,"[",1),sapply(inList,"[",2))


#' Add Column Names
#'
#' functionalizes colnames(matrix)<-newnames operation.
#' @param matix The data that the new column names are being added to
#' @param colnames The new column names
#' @return matrix with the same values but new col.name[s]
#' @examples
#' addColnames(cbind(c(1,2,3),c(4,5,6)),c("low","high"))
#' @seealso \code{\link{addRownames} \link{addNames}}
#' @export
addColnames<-function(matrix,colnames){
    colnames(matrix)<-colnames
    matrix
}

#' Add Row Names
#'
#' functionalizes rownames(matrix)<-newnames operation.
#' @param matix The data that the new row names are being added to
#' @param colnames The new row names
#' @return matrix with the same values but new row.name[s]
#' @examples
#' addRownames(cbind(c(1,2,3),c(4,5,6)),c("a","b","c"))
#' @seealso \code{\link{addColnames} \link{addNames}}
#' @export
addRownames<-function(matrix,rownames){
    rownames(matrix)<-rownames
    matrix
}

#' Add Names
#'
#' Adds row and column names to a matrix
#' @param matrix input matrix or list
#' @param colnames new colnames
#' @param rownames new rownmaes
#' @param list if not NULL matrix is a list
#' @return a list with new names
#' @examples
#' addNames(cbind(c(1,2,3),c(4,5,6)),c("high","low"),c("a","b","c"))
#' addNames(list(1,2,3),c("a","b","c"),list=TRUE)
#' @export
#' @seealso \code{\link{addColnames} \link{addRownames}}
addNames<-function(matrix,colnames,rownames=colnames,list=NULL){
    if(is.null(list)){
        colnames(matrix)<-colnames
        rownames(matrix)<-rownames
    }
    else{
        names(matrix)<-colnames
    }
    matrix        
}


#' Normalize
#'
#' privides unit normalization for numeric x
#' @param x numeric input
#' @return numeric of length x
#' @examples
#' x<-rnorm(100,10,2)
#' y<-normalize(x)
#' print(c(mean(y),sd(y)))
#' @export
normalize<-function(x){
    (x-mean(x))/(sqrt(var(x)))
}

#' Make Logical
#' 
#' transforms a vector of indicies into logical values.
#' Requires the full size of the output list
#' @param loc a list of indicies
#' @param size the length of the output vector
#' @return a vector of logicals
#' @examples
#' i<-c(1,2,3)
#' all(makeLogic(i,4) == c(TRUE,TRUE,TRUE,FALSE))
#' @export
makeLogic<-function(loc,size){
    x=rep(FALSE,size)
    x[loc]<-TRUE
    x
}


#' stem
#'
#' Produces a stem plot similar to matlab using the plotting utility.
#' @author Matti Pastell
#' @seealso \code{\code{plot}}
#' @references \url{http://www.r-bloggers.com/matlab-style-stem-plot-with-r/} 
#' @examples
#' stem(1:10,runif(10))
#' @export
stem <- function(x,y,pch=16,linecol=1,clinecol=1,...){
    if (missing(y)){
        y = x
        x = 1:length(x) }
    plot(x,y,pch=pch,...)
    for (i in 1:length(x)){
        lines(c(x[i],x[i]), c(0,y[i]),col=linecol)
    }
    lines(c(x[1]-2,x[length(x)]+2), c(0,0),col=clinecol)
}


#' ors
#'
#' Recursivly finds the union of a list of logicals
#' @seealso ~\code{\link{orM} \link{"|"}}
#' @param inList a list of vectors
#' @return the union of the list of vectors
#' @examples
#' result<-ors(list(c(1,2,3),c(3,4,5)))
#' result<-c(TRUE,TRUE,TRUE)
#' @export
ors<-function(inList){
    l<-length(inList)
    if(l==0)
        return(NULL)
    if(l==1)
        return(inList[[1]])
    if(l==2)
        return(inList[[1]] | inList[[2]])
    if(l>2){
        ltemp<-append(list(inList[[1]] | inList[[2]]),inList[3:l])
        return ( and(ltemp))
    }
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
#' @export
orM<-function(mat){
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

#' union N
#' 
#' Recursivly tales the union of a list of vectors. Takes 0 to N vectors. 
#' The default behavior for a length of 0 is to return c(), i.e. NULL.
#' unionN(a,b,c) ~ Reduce(union,list(a,b,c))
#' @param ... a variable number of vectors to have their union taken
#' @return a vector containing the union of the values.
#' @seealso \code{\link{intersectN} \link{union}}
#' @keywords union recursive
#' @examples
#' intersectN()
#' # > NULL
#' intersectN(c(1,2,3))
#' # > 1 2 3 
#' intersectN(c(1,2,3),c(2,3,4))
#' # > 2 3
#' intersectN(c(1,2,3),c(2,3,4),c(3,4,5))
#' #> 3
#' @export
unionN<-function(...){
    a<-list(...)
    len<-length(a)
    if( len==0)
        c()
    else if(len==1)
        a[[1]]
    else if (len==2)
        union(a[[1]],a[[2]])
    else if (len >2 )
        unlist(do.call(unionN,append(list(as.list(union(a[[1]],a[[2]]))),a[3:len])))
}

#' instersect N
#'
#' Recursivly intersects a list of vectors. Takes 0 to N vectors. 
#' The default behavior for a length of 0 is to return c(), i.e. NULL.
#' intersectN(a,b,c) ~ Reduce(intersect,list(a,b,c))
#' @param ... a variable number of conses to be intersected
#' @return a vector containing the intersecting values.
#' @seealso \code{\link{unionN}\link{intersect}}
#' @keywords intersect recursive
#' @examples
#' intersectN()
#' # > NULL
#' intersectN(c(1,2,3))
#' # > 1 2 3 
#' intersectN(c(1,2,3),c(2,3,4))
#' # > 2 3
#' intersectN(c(1,2,3),c(2,3,4),c(3,4,5))
#' #> 3
#' @export
intersectN<-function(...){
    a<-list(...)
    len<-length(a)
    if( len==0)
        c()
    else if(len==1)
        a[[1]]
    else if (len==2)
        intersect(a[[1]],a[[2]])
    else if (len >2 )
        unlist(do.call(intersectN,
                       append(list(as.list(intersect(a[[1]],a[[2]]))),
                              a[3:len])))
}


#' Overlap
#'
#' Takes a matrix of logical inputs, each column representing a data set each
#' row representing a location on the genome,
#'
#' @param f n by m logical matrix
#' @return A m by m matrix of column overlaps
#' @export
#' @examples
#' lmatrix<-matrix(as.logical(c(1,1,1,1,1,1,0,0,1,0,1,0)),ncol=3)
#' overlap(lmatrix)
overlap<-function(f){
    permutations<-function(lin)
        as.vector(outer(lin,lin,function(x,y) paste(x,y,sep="")))
    overlapNorm<-function(x,f){
        a<-which(f[,x[1]]==1)
        b<-which(f[,x[2]]==1)
        length(intersect(a,b))/length(union(a,b))
    }
    width<-dim(f)[2]
    m<-matrix(apply(do.call(cbind,lapply(strsplit(permutations(seq(width)),split=""),as.numeric)),2, overlapNorm,f),ncol=width)
    colnames(m)<-colnames(f)
    rownames(m)<-colnames(f)
    m
}

#' map zip list
#'
#' Takes a list of lists and groups the nth member of each list
#' into a new one.
#' @param ... an arbitraraly long set of equal length lists
#' @export
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
mapziplist<-function(...){
    inlist<-list(...)
    l<-length(inlist[[1]])
    zip<-function(i)
        unlist(lapply(inlist,"[[",i))
    lapply(seq(l),zip)
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
sortDataFrame<-function(dd,...)    
    dd[do.call(order,lapply(list(...),function(x) dd[x])),]


#' method or
#' 
#' iterates through the members of the input vector and returns the first
#' variable that has a non NULL value
#' @param ... variable length list of inputs
#' @return the first non NULL value
#' @examples
#' mor() # > NULL
#' mor(1,2,3) # > 1
#' mor(NULL,1) # > 1
#' # useful for generalizing dim(x)[2]
#' mor(dim(c(1,2,3))[2],1) #> 1
mor<-function(...){
    for(i in list(...))
        {
            if(! is.null(i))
                return(i)
        }
}

#' m Width
#'
#' Returns the width of a vector, if it is a non matrix
#' it returns 1 rather than NULL.
#' @param mat input matrix
#' @return the width of the matrix
#' @seealso \code{\link{mLength}}
#' @examples
#' mWidth(c(1,2,3)) # > 1
#' mWidth(cbind(c(1,2,3),c(4,5,6))) # > 2
mWidth<-function(mat)
    mor(dim(mat)[2],1)

#' m Length
#'
#' Determines the length of a matrix or vector.
#' 
#' @param mat input matrix
#' @return the width of the matrix
#' @seealso \code{\link{mLength}}
#' @examples
#' mLength(c(1,2,3)) # > 3
#' mLength(cbind(c(1,2,3),c(4,5,6))) # > 3
mLength<-function(mat)
    mor(dim(mat)[1],length(mat))

#' modulous
#'
#' Turns the %% infix operater into a function
#' @param x numerator 
#' @param m denominator
#' @return an numeric value representing the modlulous
#' @export
modulous<-function(x,m){
    # t1<-floor(x/m)
    # (x-t1*m)
    x %% m
}
