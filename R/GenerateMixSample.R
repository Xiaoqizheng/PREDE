#################################################################################
#Generate_bulk: generate mixed samples based on the celllines before simulation
#################################################################################

generate_bulk <- function(celltypes,nSample =100,csd = 0.1){
  k = ncol(celltypes)
  m = nrow(celltypes)
  
  ## generate composition matrix H
  H.dat = t(rdirichlet(nSample,c(rep(1,k))))
  H.dat = matrix(runif(k*nSample,min = 0,max = 1),nrow = k)
  H = apply(H.dat,2, function(x) x/sum(x))
  Y = as.matrix(celltypes) %*% H
  
  ## add noise
  sd = as.vector(Y)*csd 
  noise = matrix(rnorm(n = m*nSample,mean = 0,sd = sd),nrow = m)
  Y = Y + noise
  
  Y[Y<0] = 0
  
  return(list(Y = Y,H = H,W = celltypes))
}





