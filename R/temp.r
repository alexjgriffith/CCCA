# this file is for temporary c-wrapping
# it is dependant on the inclusion of mulcal

#' @export
regionWrapper<-function(temp,chroms){
    width<-length(chroms)
    buffer<-as.ineger(rep(0,width))
    outMatrix<-integer((width+2)*(length-1))
    length<-length(temp[,1])
    .C("region",temp[,1],temp[,3],temp[,2],buffer,outMatrix,length,width)
}

#' @export
peakRegions<-function(data){
    chroms<-levels(data[,1])    
    subset<-data[,data[,1]==chroms[1]]
    temp<-rbind(cbind(subset[,2],subset[,"name"],0),
                cbind(subset[,3],subset[,"name"],1))
    regionWrapper(temp,chroms)
}
