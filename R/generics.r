#' Contribution Plot
contribution<-function(x,...){
    UseMethod("contribution",x)
}

#' Contribution Plot
clust<-function(x,...){
    UseMethod("clust",x)
}

#' Normalize Plot
normalize<-function(x,...){
    UseMethod("normalize",x)
}

#' AddReg Plot
addReg<-function(x,...){
    UseMethod("addReg",x)
}

#' AddReg Plot
addFasta<-function(x,...){
    UseMethod("addFasta",x)
}


#' @rdname contribution
#' @method contribution default
#' @export
contribution.default<-function(x,...){
    stop("Currently only implemented for prc class.")
}

#' @rdname clust
#' @method clust default
#' @export
clust.default<-function(x,...){
    hclust(x,...)
}

#' @rdname normalize
#' @method normalize default
#' @export
normalize.default<-function(x,...){
    
}

#' @rdname addReg
#' @method addReg default
#' @export
addReg.default<-function(){
}

#' @rdname addFasta
#' @method addFasta default
#' @export
addFasta.default<-function(){
}
