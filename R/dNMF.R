#' discrimant Non-Negative Matrix Factorization (dNMF)
#' 
#'  The dNMF algorithm with the additional Fisher criterion on the cost 
#'  function of conventional NMF was designed to increase class-related
#'  discriminating power.

dat = rbind(matrix(c(rep(3, 16), rep(8, 24)), ncol=5), 
            matrix(c(rep(5, 16), rep(5, 24)), ncol=5),
            matrix(c(rep(18, 16), rep(7, 24)), ncol=5)) + 
    matrix(runif(120,-1,1), ncol=5)

trainlabel <- c(1,1,2,2,2)
res <- dNMF(dat, trainlabel, r=2, lambada = 1)
res$H

dNMF <- function(dat, trainlabel, r=2, lambada=0.1, maxIter=1000, tol=1e-7, verbose=TRUE){
    
    dat <- as.matrix(dat)
    nFea = nrow(dat); nSmp = ncol(dat)
    eps = .Machine$double.eps
    
    #init the H0 and W0 matrix
    H = matrix(runif(r*nSmp, eps), r, nSmp)
    H = pmax(H,0)
    W = matrix(runif(nFea*r, eps), nFea, r)
    W = W/colSums(W)
    W = pmax(W, eps)
    b = pmax(abs(W %*% H), eps)
    
    #calculate EL distance of two matrix
    for(i in unique(trainlabel)){
        Hclass[,which(trainlabel==i)] = matrix( rep(rowMeans(H[,which(trainlabel==i)]), length(which(trainlabel==i))), r, length(which(trainlabel==i)))
    }
    Hsum = matrix(rep(rowMeans(H), nSmp),r, nSmp)
    Nk = matrix(rep(table(trainlabel), times=r * table(trainlabel)), nrow=r, byrow=FALSE)
    Ed = lambada * sum(Nk * (Hclass -Hsum)^2) / (2 * nSmp * r)
    
    obj0 = norm(dat-b+eps, type="F") - Ed
    obj1 =obj0
    
    #averaing matrix Ma over all samples and Mc for wach class
    Ma = matrix(1/nSmp, nrow=nSmp, ncol=nSmp)
    Mc <- list()
    for (i in unique(trainlabel)){
        Mc[[i]] = matrix(1/table(trainlabel)[i], table(trainlabel)[i], table(trainlabel)[i])
    }
    Mc = Matrix::bdiag(Mc)
    
    Hclass = H
    final = Inf
    count = 1
    obj_stack =  vector (mode="numeric", length = maxIter)
    
    while (final > tol && count <= maxIter) {
        
        # update W and H
        W = W * (dat %*% t(H) / (W %*% (H %*% t(H))) )
        lambadaX = lambada * nFea/r * H
        H = H * ( (t(W) %*% dat + lambadaX %*% Mc) / (t(W) %*% W %*% H + lambadaX %*% Ma) )
        
        # H normalization (Formula 11)
        H = as.matrix(H)
        H = H / (rowSums(H)/nSmp + eps)
        
        # convergence
        for(i in unique(trainlabel)){
            Hclass[,which(trainlabel==i)] = matrix( rep(rowMeans(H[,which(trainlabel==i)]), length(which(trainlabel==i))), r, length(which(trainlabel==i)))
        }
        Hsum = matrix(rep(rowMeans(H), nSmp),r, nSmp)
        Nk = matrix(rep(table(trainlabel), times=r * table(trainlabel)), nrow=r, byrow=FALSE)
        Ed = lambada * sum(Nk * (Hclass -Hsum)^2) / (2 * nSmp * r)
        
        b <- pmax(abs( W%*%H ), eps)
        obj2 = obj1
        obj1 = norm(dat-b+eps, type="F") - Ed
        final = abs(obj1-obj2) / (abs(obj1-obj0) + eps)
        
        obj_stack[count] = obj1
        count = count + 1
        
    }
    list(V=dat, W=W, H=H, trainlabel=trainlabel, count=count,
         final=final, obj_stack=obj_stack, r=r, call=match.call())
    
}


