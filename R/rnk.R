#' write rnk to a file from matrix W.
#' 
#' write a rnk file from matrix W in a returned object of function \code{DNMF}.
#' The rnk format is referred \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{RNK}
#' 
#' @param object a returned object of function \code{DNMF}
#' @param fn the output filename. Default is "./tmp.rnk"
#' @export
#' @examples
#' \dontrun{
#' rnk(dnmf_result, fn="./tmp.rnk")
#' }
rnk <- function (object, fn="./tmp.rnk"){
    write.table(object$rnk, fn, sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
}
