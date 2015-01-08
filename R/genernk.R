#' ranking the genes via DNMF
#' 
#' ranking the genes via DNMF from a returned object of function \code{DNMF}.
#' 
#' @param object, a returned object of function \code{DNMF}
#' @export
#' @examples
#' \dontrun{
#' genernk(dnmf_result)
#' }
genernk <- function (object){
    # To distinguish the row of H, assume the 1st row of H is down_regulated (relative the Ctrl), if not, make it to be.
    if (object$H[1,1] < object$H[1,ncol(object$H)]) {
        tmp_W = object$W[,2]
        object$W[,2] = object$W[,1]
        object$W[,1] = tmp_W
        rm (tmp_W)}
    return (object$W[,2]-object$W[,1])
    
}