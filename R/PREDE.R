#################################################################################
# PREDEModel: Partial-reference based deconvolution model
#################################################################################
## PREDE: partial reference-based deconvolution model
## Y: the profiles of tumor samples, they can be gene expression proflies,methylation profile, etc.
## W: the profiles of cell types
## W1: the profiles of partial reference cell lines
## type: "GE" for gene; "ME" for methylation
## K: the total number of cell types including both known and unknown cell types

PREDE <- function(Y,W = NULL,W1,type = "GE",K,iters = 500,rssDiffStop=1e-10){
    # Do a linear transformation of Y to fit the QPfunction function
    m = max(Y)
    Y = Y*100/m
    if (is.null(W1)){
    ### ref-free
    print("This is ref-free case")
    out.RF = RF(Y,W = W, K = K,type = type)
    W.pred = (out.RF$W)*m/100
    H.pred = out.RF$H

    if(is.null(W)){
      ## output result directly
      out = list(W=W.pred, H = H.pred)
    } else {
      ## Adjust W.pred and H.pred using W
      out = adjustWH(W,W.pred,H.pred)
    }
    return(out)
  } else {
    W1=as.matrix(W1)
    W1 = W1*100/m
    if (ncol(W1) == K){
      #### ref-based
      print("This is ref-based case")
      return(RB(Y,W1))
    } else if (ncol(W1) < K){
      ### partial ref
      mu02 <- RefFreeCellMixInitialize(Y,K=K-ncol(W1),method="ward.D2")
      mu02 <- as.matrix(mu02)
      mu0 <-cbind(as.matrix(W1),mu02)
      rss0 = 0

      for(i in 1:iters){
        flag <- !apply(is.na(mu0),1,any)
        omega <- QPfunction(Y[flag,],mu0[flag,])
        omega1 <- omega[,1:ncol(W1)]
        omega2 <- omega[,(ncol(W1)+1):ncol(mu0)]

        if (type == "GE"){
          mu1 <- QPfunction(t(Y-W1%*%t(omega1)), omega2, sumLessThanOne=FALSE)
        } else if (type == "ME"){
          mu1 <- QPfunction(t(Y-W1%*%t(omega1)), omega2, sumLessThanOne=FALSE,W_lessThanOne = T)
        } else {
          stop("Please input the correct data type!")
        }
        mu = cbind(W1,mu1)
        rss.new = norm(Y - mu%*%t(omega),"F")^2
        ## check convergence
        dd = abs(rss.new-rss0)
        if(dd < rssDiffStop)
          break

        ## updata
        mu0 = cbind(W1,mu1); rss0 = rss.new
      }

      W.pred = mu*m/100
      H.pred = t(omega)

      if(is.null(W)){
        out = list(W=W.pred, H = H.pred)
      } else {
        # Adjust W.pred and H.pred using W
        out = adjustWH(W,W.pred,H.pred)
      }

      return(out)
    } else {
      stop("ncol(W1) should be less than or equal to ncol(W)!")
    }
  }
}


RB <- function(Y,W){
  H.pred = t(QPfunction(Y = Y, Xmat = W, sumLessThanOne=TRUE, nonNeg=!sumLessThanOne))
  return(list(H = H.pred))
}

RF <- function(Y,W = NULL,K,type = "GE",iters = 500,rssDiffStop=1e-10){
  mu0 <- RefFreeCellMixInitialize(Y,K=K,method="ward.D2") ## initial profiles of cell types
  mu0 <- matrix(mu0,ncol = K)
  rss0 = 0
  for(i in 1:iters){
    flag <- !apply(is.na(mu0),1,any)
    omega <- QPfunction(Y[flag,],mu0[flag,],sumLessThanOne=TRUE)
    if (type == "GE"){
      mu <- QPfunction(t(Y), omega, sumLessThanOne=FALSE,W_lessThanOne = F)
    } else if (type == "ME"){
      mu <- QPfunction(t(Y), omega, sumLessThanOne=FALSE,W_lessThanOne = T)
    } else {
      stop("Please input the correct data type!")
    }

    rss.new = norm(Y - mu%*%t(omega),"F")^2
    ## check convergence
    dd = abs(rss.new-rss0)
    if(dd < rssDiffStop)
      break

    ## updata
    mu0 = mu; rss0 = rss.new
  }

  W.pred = mu
  H.pred = t(omega)

  if(is.null(W)){
    out = list(W=W.pred, H = H.pred)
  } else {
    # Adjust W.pred and H.pred using W
    out = adjustWH(W,W.pred,H.pred)
  }
  out
}

QPfunction <- function(Y, Xmat, sumLessThanOne=TRUE, nonNeg=!sumLessThanOne, W_lessThanOne = F){
  Xmat = as.matrix(Xmat)
  nCol = dim(Xmat)[2]   # the number of cell types
  nSubj = dim(Y)[2]     # sample size
  mixCoef = matrix(0, nSubj, nCol) ## the proportion of cell types in each samples
  rownames(mixCoef) = colnames(Y)
  colnames(mixCoef) = colnames(Xmat)

  # calculate H matrix from W
  if(sumLessThanOne){
    Amat = cbind(rep(-1,nCol), diag(nCol))
    b0vec = c(-1,rep(0,nCol))

    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i]))
      Dmat = t(Xmat[obs,]) %*% Xmat[obs,]
      dvec = t(Xmat[obs,])%*%Y[obs,i]
      meq = 1
      mixCoef[i,] = solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=b0vec, meq=meq)$sol
    }
  }

  # calculate W matrix from H
  else if(nonNeg){

    if(W_lessThanOne == F){ ## only x >= 0, for gene expression
      Amat = cbind(diag(nCol))
      b0vec = c(rep(0,nCol))
    } else if(W_lessThanOne == T){ ## 0 <= x <= 1, for DNA methylation
      Amat = cbind(-diag(nCol), diag(nCol))
      b0vec = c(rep(-1,nCol),rep(0,nCol))
    }

    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i]))
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec)$sol
    }
  }

  return(mixCoef)
}


adjustWH <- function(W,W.pred,H.pred){
  MAT = cor(W,W.pred)
  n = min(dim(MAT))
  a.index = 1:n

  for(i in 1:n) {
    z = which(MAT == max(MAT),arr.ind=T)
    a.index[z[1,1]] = z[1,2]
    MAT[,z[1,2]] = -1
    MAT[z[1,1],] = -1
  }

  W.adj = W.pred[,a.index]
  H.adj = H.pred[a.index,]
  return(list(W = W.adj,H = H.adj))
}

RefFreeCellMixInitialize <- function(Y,K=2,Y.Distance=NULL, Y.Cluster=NULL,
                                     largeOK=FALSE, dist.method = "euclidean", ...){

  if(!is.matrix(Y) | !is.numeric(Y)){
    stop("Y is not a numeric matrix\n")
  }
  n <- dim(Y)[2]

  if(is.null(Y.Cluster)){
    if(is.null(Y.Distance)){
      if(n>2500 & !largeOK){
        stop("Y has a large number of subjects!  If this is what you really want, change 'largeOK' to TRUE\n")
      }
      Y.Distance <- dist(t(Y),method=dist.method)
    }
    Y.Cluster <- hclust(Y.Distance,...)
  }

  classes <- cutree(Y.Cluster, K)
  s <- split(1:n,classes)

  sapply(s, function(u) apply(Y[,u,drop=FALSE],1,mean,na.rm=TRUE))
}




