#' Discrimant Non-negative Matrix Factorization.
#'
#' The DNMF is to extend the Non-negative Matrix Factorization algorithm in 
#' order to extract features that enforce not only the spatial locality, but
#'  also the separability between classes in a discriminant manner.
#' 
#' The main algorithm is based on 
#' \href{http://www.ncbi.nlm.nih.gov/pubmed/16722172}{Zafeiriou, S., et al. 
#' (2006) Exploiting discriminant information in 
#' nonnegative matrix factorization with application to frontal face 
#' verification, IEEE transactions on neural networks, 17, 683-695}, 
#' with some \strong{CORRECTIONs}. 
#'
#' @param data a matrix, like expression profilings of some samples. the columns are samples and the rows are gene's expression.
#' @param trainlabel a numeric vector of sample type of all the samples, this vector should ONLY contain 1 and 2 so far and length of it should equal the column (sample) size of data.
#' @param r the dimension of expected reduction dimension, with the default value 2.
#' @param gamma the tradeoff value for the within scatter matrix, with the default value 0.1.  
#' @param delta the tradeoff value for the between scatter matrix, with the default value 1e-4.
#' @param maxIter the maximum iteration of update rules, with the default value 1000.
#' @param tol the toleration of coverange, with the default value 1e-7.
#' @param plotit whether plot H (V=WH). Default: FALSE, other options: TRUE.
#' @export
#' @examples
#' r =2
#' data =  matrix(1:12,3,4)
#' trainlabel = c(1,1,2,2)
#' DNMF_result <- DNMF(data, trainlabel, r)

DNMF <- function(data,trainlabel, r=2, gamma=0.1, delta=0.0001, maxIter=1000, 
                 tol=1e-7, plotit=FALSE) {
	
    data <- as.matrix(data)
    nFea = nrow(data); nSmp = ncol(data)
    eps = .Machine$double.eps
    
    #init the H0 and W0 matrix
    H = matrix(runif(r*nSmp, eps), r, nSmp)
    H = pmax(H,0)
    W = matrix(runif(nFea*r, eps), nFea, r)
    W = W/colSums(W)

    #calculate KL divergence of two matrix
    b = pmax(abs(W %*% H), eps)
    obj0 = sum(data*log((data+eps)/(b-eps))-data+b)
    obj1 = obj0 
    
    ##########################
    E = matrix(1, nFea, nSmp)
    # N is just 1/Nr in paper, the weighted matrix
    SmpCount_withinClass = vector(mode="numeric", length(trainlabel))
    SmpCount_withinClass = as.vector(table(trainlabel)[trainlabel])
    N = matrix(1/SmpCount_withinClass, r, nSmp, byrow=T)     

    final = Inf
    count = 1
    Hclass = H  
    obj_stack =  vector (mode="numeric", length = maxIter)
    
while (final > tol && count <= maxIter) {

    #update H with the objective function includes KL divergence

    for(i in unique(trainlabel)){
        Hclass[,which(trainlabel==i)] = matrix( rep(rowSums(H[,which(trainlabel==i)]), length(which(trainlabel==i))), r, length(which(trainlabel==i)))
    }
    Hclass = Hclass - H
    Hsum = matrix(rep(rowSums(H),ncol(H)),r, ncol(H)) - H

    tmp_a = 4*gamma + 4*(delta/nSmp - (gamma+delta)*N)   #2a
    tmp_b = 1 + 2*delta/nSmp * Hsum - 2*(delta+gamma)*N * Hclass
    tmp_c = - t(W) %*% (data / ( W %*% H + eps) ) * H
    H = (sqrt(tmp_b^2 - 2*tmp_a*tmp_c) - tmp_b ) /tmp_a
    H = pmax(H, eps)  

    #######################################
    #update W
    W = (W/(E%*%t(H)))*(data/(W%*%H)) %*% t(H)
    H = diag(colSums(W)) %*% H
    W = W/ matrix(rep(colSums(W), each=nrow(W)), nrow(W), ncol(W))
    W = pmax(W,eps)
    
    obj2 = obj1
    b = pmax(abs(W%*%H),eps)
    obj1 = sum( data*log((data+eps)/(b-eps)) - data + b )
    final = abs(obj1-obj2) / abs(obj1-obj0)
    #obj_stack[count] = final
    obj_stack[count] = obj1
    count = count + 1
}
    # to plot the convergence of the object function
    if (plotit){
    	heatmap(H)
    }
    

    # check H. Setting down-regulted metagene in row 1 of H, 
    # while up-regulated metagene in row 2 of H.
    l1 <- apply(H[,which(trainlabel == names(table(trainlabel))[1])], 1, mean)
    l2 <- apply(H[,which(trainlabel == names(table(trainlabel))[2])], 1, mean)
    meanH <- cbind(l1, l2)
    if (l1[1] < l2[1] && l1[2] > l2[2]) {
        H <- rbind(H[2,], H[1,])
        W1 <- cbind(W[,2], W[,1]) 
    } else if ((l1[1] < l2[1] && l1[2] < l2[2]) || (l1[1] > l2[1] && l1[2] > l2[2])) {
        stop("Failed. Run DNMF again after restart R.")
    }
    print ("The label-specific mean of H:")
    print (meanH)

    list(W=W, H=H, trainlabel=trainlabel, meanH=meanH, delta=delta, gamma=gamma, count=count,
         final=final, obj_stack=obj_stack, r=r, call=match.call())
}


