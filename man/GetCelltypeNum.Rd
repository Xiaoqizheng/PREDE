\name{GetCelltypeNum }
\alias{GetCelltypeNum }
\title{
get optimal number of total cell types in the mixture tumor sample
}
\description{
get the optimal number of total cell types by computing AIC value.
}
\usage{
GetCelltypeNum(Y,W = NULL,W1 = NULL,maxK = 50)
}
\arguments{
  \item{Y}{matrix of profiles of mixture samples, which can be gene expression proflies,methylation        profile, etc. }
  \item{W}{matrix of the profiles of cell types.}
  \item{W1}{matrix/vector of the profiles of partial reference cell types.}
  \item{maxK}{an upper bound of the number of total cell types. }
}
\details{
Before the deconvolution of mixture sample into distinct cell populations based on the partial reference, the number of total cell types should be specified. We obtain the optimal number of total cell types by computing the Akaike information criterion (AIC). Given the range of K, we get a series of values of AIC and take the number of total cell types with minimum aic as the optimal number of total cell types.
}
\value{
a vector for AIC values under different numbers of total cell types, and a value for optimal number of total cell types .
}
\references{
Y. Qin, W. Zhang, S. Nan, N. Wei and X. Zheng (2020). Deconvolution of heterogeneous tumor samples using partial reference signals. Submitted.
}
\author{
Xiaoqi Zheng \email{xqzheng@shnu.edu.cn}.
}
\examples{

## load example data
data(lung_exp)
W <- lung_exp[,1:6]

## generate the mixed samples
bulk <- generate_bulk(W,nSample =100,csd = 0.1)

## select the feature
feat <- select_feature(mat = bulk$Y,method = "cv",nmarker = 1000,startn = 0)

## get optimal number of total cell types
OptimalK <- GetCelltypeNum(bulk$Y[feat,],W=NULL,W1=W[feat,1:4],maxK = 10)
}
