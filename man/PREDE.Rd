\name{PREDE}
\alias{PREDE}
\title{
Partial-reference based deconvolution model
}
\description{
Estimate the profiles of all cell types and the proportion of cell types in tumor samples
}
\usage{
PREDE(Y,W = NULL,W1,type = "GE",K,iters = 500,rssDiffStop=1e-10)
}
\arguments{
  \item{Y}{matrix of profiles of tumor samples, which can be gene expression proflies,methylation profile, etc. }
  \item{W}{matrix of the profiles of cell types.}
  \item{W1}{matrix/vector of the profiles of partial reference cell types.}
  \item{type}{type of tumor sample data. Options are "GE","ME", etc.}
  \item{K}{number of cell types.}
  \item{iters}{number of iterations to execute.}
  \item{rssDiffStop}{conditions for iteration termination.}
}
\details{
Partial-reference based decomposition model could estimate the profiles of unknown cell types and the proportions of cell types. The number of total cell types K should be specified when running the function PREDE. The optimal value of K could be determined by the output value AIC of the function PREDE. If the profile of cell types W is non-empty, we can use adjustWH to adjust the order of predicted cell types.
}
\value{
a matrix for profiles of cell types,a matrix for the proportions of cell types and a value of AIC.
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

## Partial-reference based deconvolution without total cell types W
bulk <- GenerateMixSample(W,nSample =100,csd = 0.1)
PREDE(bulk$Y,W1=W[,1:4],type = "GE",K=7,iters = 500,rssDiffStop=1e-10)


## Partial-reference based deconvolution with total cell types W
bulk <- GenerateMixSample(W,nSample =100,csd = 0.1)
PREDE(bulk$Y,W,W1=W[,1:4],type = "GE",K=7,iters = 500,rssDiffStop=1e-10)

}
