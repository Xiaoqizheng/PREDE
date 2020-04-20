
#################################################################################
## select optimal number of total cell types by computing AIC value
#################################################################################
GetCelltypeNum <- function(Y,W = NULL,W1 = NULL,maxK = 50){
    W1 = as.matrix(W1)
    AIC = c()
    for (k in (ncol(W1)+1):maxK){
        PR = PREDE(Y,W = W,W1 = W1,type = "GE",k,iters = 100,rssDiffStop=1e-4)
        W.pred = PR$W
        H.pred = PR$H
        rss = norm(Y - W.pred %*% H.pred,type = "F")^2
        nSamples = ncol(Y)*nrow(Y) ### sample size * number of gene
        nParam = k*(nrow(Y)+ncol(Y)) - nrow(W1)*ncol(W1) ### total number of parameters
        aic = nSamples*log(rss/nSamples)+ 2*nParam + (2*nParam*(nParam+1))/(nSamples-nParam-1)
        AIC = rbind(AIC,aic)
    }
    ##select the optimal number of cell types
    K = which.min(AIC)+ncol(W1)
    out = list(AIC=AIC,K=K)
    out$AIC = AIC
    out$K = K
    return(out)
}


