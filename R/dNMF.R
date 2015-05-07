
dat = rbind(matrix (rep(c(3,8, 12),each=16),8), 
            matrix (rep(c(5,5,9),each=8),4),
            matrix (rep(c(18,7, 12),each=16),8)) + 
    matrix(runif(120,-1,1), ncol=6)

dat = rbind(matrix(c(rep(3, 16), rep(8, 24)), ncol=5), 
            matrix(c(rep(5, 16), rep(5, 24)), ncol=5),
            matrix(c(rep(18, 16), rep(7, 24)), ncol=5)) + 
    matrix(runif(120,-1,1), ncol=5)

trainlabel <- c(1,1,2,2,2)

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
    obj0 = norm(dat-b+eps, type="F")
    obj1 =obj0
    
    #averaing matrix Ma over all samples and Mc for wach class
    Ma = matrix(1/nSmp, nrow=nSmp, ncol=nSmp)
    Mc <- list()
    for (i in unique(trainlabel)){
        Mc[[i]] = matrix(1/table(trainlabel)[i], table(trainlabel)[i], table(trainlabel)[i])
    }
    Mc = Matrix::bdiag(Mc)
    
    final = Inf
    count = 1
    obj_stack =  vector (mode="numeric", length = maxIter)
    
    while (final > tol && count <= maxIter) {
        W = W * (dat %*% t(H) / (W %*% H %*% t(H)) )
        lambada = lambada * nFea/r * H
        H = H * ( (t(W) %*% dat + lambada %*% Mc) / (t(W) %*% W %*% H + lambada %*% Ma) )
        
        # H normalization (Formula 11)
        H = as.matrix(H)
        H = H / (rowSums(H)/nSmp)
        
        b <- pmax(abs( W%*%H ), eps)
        obj2 = obj1
        obj1 = norm(dat-b+eps, type="F")
        final = abs(obj1-obj2) / abs(obj1-obj0)
        
        obj_stack[count] = obj1
        count = count + 1
    }
    list(V=dat, W=W, H=H)
    
}

res <- dNMF(dat, trainlabel, r=2, lambada = 0.01)
res$H

