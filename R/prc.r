#' PCA
#'
#' Normalizes a data matrix and then preformes quantile normalization
#' using \code{prcomp} on the transform of the normalized data. The
#' eigenvectors are returned in a list with the normalized data
#' @param data The matrix to be analyzed
#' @param norm The normalization method
#' norm may either be a user provided function or one of the
#' following
#' \describe{
#' \item{rowSumOne}{ all rows sum to 1}
#' \item{colSumOne}{ all columns sum to 1}
#' \item{normRow}{all rows have mean 0 and sd 1}
#' \item{normCol}{all cols have mean 0 and sd 1}
#' \item{rowsVarOne}{ all rows have sd 1}
#' \item{colsVarOne}{ all cols have sd 1}
#' \item{qn}{ quantile normalization}
#' \item{none}{ no normalization}
#' }
#' @return list($normData,$eigenVectors)
#' @export
makrPRC<-function(data,norm="qn"){
    #print(head(data))
    if( is.function(norm))
       normData<-norm(data)
    else
        normData<-switch(norm,
               rowSumOne=t(apply(data,1, function(x) {x/sum(x)})),
               colSumOne=apply(data,2, function(x) x/sum(x)),
               normRow=t(apply(data,1, function(x) (x-mean(x))/var(x))),
               normCol=apply(data,2, function(x) (x-mean(x))/var(x)),
               rowsVarOne=t(apply(data,1, function(x) x/var(x))),
               colsVarOne=apply(data,2, function(x) x/var(x)),
               qn=CCCA::._qn(data),
               none=data,
               data)
    #print(head(normData))
    prc<-prcomp(t(normData))$rotation
    list(normData=normData,eigenVectors=prc)
}
