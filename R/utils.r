or<-function(...){
    for(i in list(...))
        {
            if(! is.null(i))
                return(i)
        }
}

mWidth<-function(x)
    or(dim(x)[2],1)

mLength<-function(x)
    or(dim(x)[1],length(x))
    
and<-function(inList){
    l<-length(inList)
    if(l==0)
        return(NULL)
    if(l==1)
        return(inList[[1]])
    if(l==2)
        return(inList[[1]] & inList[[2]])
    if(l>2){
        ltemp<-append(list(inList[[1]] & inList[[2]]),inList[3:l])
        return ( and(ltemp))
    }
}
