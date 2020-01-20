\name{GenerateMixSample}
\alias{GenerateMixSample}
\title{
generate the mixed samples
}
\description{
generate the mixed samples based on the cell types and the randomly generated proportion.
}
\usage{
GenerateMixSample(celltypes,nSample =100,csd = 0.1)
}
\arguments{
  \item{celltypes}{ gene expression profiles or methylation profiles of cell type.}
  \item{nSample}{the number of generated mixed samples.}
  \item{csd}{different levels of noise.}
}
\details{ This function is mostly used in simulation studies. The gene expression profiles of some cancer cell lines were selected as basis matrix, the proportion matrix H is randomly generated under the uniform distribution. The gene expression matrix of mixture samples can be obtained by multiplying W and H matrices, followed by an additional error matrix of normal distribution with mean equal to zero and different levels of noise.
}
\value{
the profile matrix of the mixed samples, the profile matrix of the cell types and the proportion matrix.
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

## generate the mixed samples based on the profile matrix of the cell types W
bulk <- GenerateMixSample(W,nSample =100,csd = 0.1)

}
