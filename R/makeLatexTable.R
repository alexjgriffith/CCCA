#' make table
#'
#' makes a latex table from a r matrix
#' @export
#' @examples
#' a<-matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),
#'           ncol=9,dimnames=list(list("alex",2),
#'           list("a","b","c","d","e","f","g","h","i") ))
#' # > a
#' #      a b c d  e  f  g  h  i
#' # alex 1 3 5 7  9 11 13 15 17
#' # 2    2 4 6 8 10 12 14 16 18
#' cat(makeTable()(a,"llllllll",label="tab:centers",caption="data table"))
#' 
#' defaultFun<-function(data,layout,cn=" ",caption=NULL,label=NULL){
#'   rownames<-rownames(data)
#'   colnames<-colnames(data)
#'   ins<-Map(function(...) unlist(list(...)),rownames,apply(data,1,list))
#'   wf<-.p(tags,"table",opts="h")
#'   if(! is.null(caption))
#'     wf<-.c(wf,.p(singleTag,"caption",caption))
#'   if(! is.null(label))
#'     wf<-.c(wf,.p(singleTag,"label",label))
#'   wf<-.c(wf,
#'   .p(tags,'center'),
#'   .p(tags,"tabular",opts=layout,n=c("{","}")))    
#'   wf(paste(makeLines(ins,paste(vb(),
#'                      vb(makeLine(unlist(list(cn,colnames)))))),vb()))
#' makeT<-makeTable(defaultFun)
#' cat(makeT(a,"llllllll",label="tab:centers",caption="data table"))
makeTable<-function(fun=defaultFun){
    # Paritialy apply a function
    # .p :: Function ...,a ->  Function a
    .p<-function(fun,...){    
        y<-list(...)
        function(...)do.call(fun,append(y,list(...)))}
    
    # combine n functions
    # .c :: Function a => a -> a
    .c<-function(...){
        ._c2<-function(a,b){function(x) a(b(x))}
        ._cAux<-function(ret,funs){
            if(length(funs)==1)
                ._c2(ret,funs[[1]])
            else
                ._cAux(.c2(ret,funs[[1]]),funs[2:length(funs)])
        }   
        funs<-list(...)    
        dess<-as.character(unlist(
            cut(length(funs),c(-1,0,1,2,Inf),
                c("none","one","two","twoOrMore"))))
        switch(dess,
               twoOrMore= ._cAux(._c2(funs[[1]],funs[[2]]),
                   funs[3:length(funs)]),
               two=._c2(funs[[1]],funs[[2]]),
               one=funs[[1]],
               none=function(x) x)           
    }
    
    makeLines<-function(manyVals,x){
        prep<-function(x)as.character(unlist(x))
        if(length(manyVals)==1)
            makeLine(prep(manyVals[1]),x)
        else
            makeLines(manyVals[2:length(manyVals)],
                      makeLine(prep(manyVals[1]),x))
    }
    
    makeLine<-function(vals,x="")
        paste(x,paste(vals,collapse=" & ")," \\\\\n",sep="")
    
    vb<-function(x="")
        paste(x,"\\hline\n",sep="")
    
    tags<-function(tag,x="",opts=NULL,n=c("[","]")){
        if(is.null(opts))
            line<-paste("\\begin{",tag, "}\n",sep="")
        else
            line<-paste("\\begin{",tag, "}",n[1],paste(opts,sep=" "),n[2],"\n",sep="")
        paste(line,x,"\\end{",tag,"}\n",sep="")
    }
    
    singleTag<-function(tag,value,x=""){
        paste(x,"\\",tag,"{",value,"}\n",sep="")
    }
    defaultFun<-function(data,layout,cn=" ",caption=NULL,label=NULL){
        rownames<-rownames(data)
        colnames<-colnames(data)
        ins<-Map(function(...) unlist(list(...)),rownames,apply(data,1,list))
        wf<-.p(tags,"table",opts="h")
        if(! is.null(caption))
            wf<-.c(wf,.p(singleTag,"caption",caption))
        if(! is.null(label))
            wf<-.c(wf,.p(singleTag,"label",label))
        wf<-.c(wf,
               .p(tags,'center'),
               .p(tags,"tabular",opts=layout,n=c("{","}")))    
        wf(paste(makeLines(ins,paste(vb(),vb(makeLine(unlist(list(cn,colnames)))))),vb()))
    }
    fun
}
