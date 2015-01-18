#' write a rnk file from matrix W.
#' 
#' write a rnk file from matrix W in a returned object of function \code{DNMF}.
#' The rnk format is referred \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RNK:_Ranked_list_file_format_.28.2A.rnk.29}{RNK}
#' 
#' @param object, a returned object of function \code{DNMF}
#' @param fn, the output filename
#' @export
#' @examples
#' \dontrun{
#' rnk(dnmf_result, fn="./tmp.rnk")
#' }
rnk <- function (object, fn="./tmp.rnk"){
    rnk <- object$W[,2] - object$W[,1]
    write.table(rnk, fn, sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
}