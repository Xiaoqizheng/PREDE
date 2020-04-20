\name{generate_bulk}
\alias{generate_bulk}
\title{
generate the mixture samples
}
\description{
generate the mixture samples based on the profile matrix of cell types and the proportion of cell types.
}
\usage{
generate_bulk(celltypes,nSample =100,csd = 0.1)
}
\arguments{
  \item{celltypes}{ gene expression profiles or methylation profiles of cell type.}
  \item{nSample}{the number of generated mixture samples.}
  \item{csd}{different levels of noise.}
}
\details{ In simulation studies, some gene expression profiles of cancer cell lines were selected as basis matrix, the proportion matrix H is randomly generated under the Dirichlet distribution. The gene expression matrix of mixture samples can be obtained by multiplying W and H matrices, followed by an additional error matrix of normal distribution with mean equal to zero and different levels of noise.
}
\value{
The function generate_bulk returns a list containing the following components:
\item{Y}{the profile matrix of the mixture samples.}
\item{W}{the profile matrix of the cell types.}
\item{H}{the proportion matrix of the cell types.}
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
bulk <- generate_bulk(W,nSample =100,csd = 0.1)

}
