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
  \item{Y}{matrix of profiles of tumor samples, which can be gene expression proflies,methylation        profile, etc. }
  \item{W}{the profile matrix of cell types.}
  \item{W1}{matrix/vector of the profiles of partial reference cell types.}
  \item{type}{type of tumor sample data. Options are "GE","ME", etc.}
  \item{K}{the number of total cell types.}
  \item{iters}{the maximum number of iterations to execute.}
  \item{rssDiffStop}{the threshold of the interations terminate. The iterations end when either the maximum number of gnerations has been reached or the absolute error is less than the given threshold 'rssDiffStop'.}
}
\details{
Partial-reference based decomposition model could estimate the profiles of unknown cell types and the proportions of cell types. The number of total cell types K should be specified when running the function PREDE. The optimal value of K could be determined by the output value AIC of the function PREDE. If the profile of cell types W is non-empty, we can use adjustWH to adjust the order of predicted cell types.
}
\value{
a matrix for profiles of cell types and a matrix for the proportions of cell types.
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
bulk <- generate_bulk(W,nSample =100,csd = 0.1)
feat <- select_feature(mat = bulk$Y,method = "cv",nmarker = 1000,startn = 0)
PREDE(bulk$Y[feat,],W1=W[feat,1:4],type = "GE",K=7,iters = 100,rssDiffStop=1e-5)

}
