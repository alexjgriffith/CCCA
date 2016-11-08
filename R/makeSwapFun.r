#' string to swap fun
#'
#' takes a string and produces a swap fun from it. The string is
#' broken allong spaces.
#' @param x input string
#' @return a function mapping the even values of x to the odd
#' @examples
#' makeSwapFun("a b c d")("b")=="a"
#' @export
makeSwapFun<-function(x){
    do.call(CCCA:::._genSwapFun,CCCA:::._splitZip(CCCA:::._createZip(strsplit(x," ")[[1]])))
}
