\name{SelectFeature}
\alias{SelectFeature}
\title{
Select feature for profile matrix of mixed samples
}
\description{
Select feature for profile matrix of mixed samples before applying PREDE model
}
\usage{
SelectFeature(mat,method,nmarker, startn=0)
}
\arguments{
  \item{mat}{matrix of profiles of mixed samples. }
  \item{method}{method of selecting features. Options are "random","cv" and "topvar". }
  \item{nmarker}{number of selecting features.}
  \item{startn}{startn+1 represents the position where the first feature starts.}
}
\details{  Since the total number of genes or CpG sites is huge compared with the number of samples and cell populations, it is necessary to select features. nmarker represents the number of selected features. startn+1 represents the position where the first feature starts.
}
\value{
The row in which the selected feature are.
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
feat = SelectFeature(mat = bulk$Y,method = "cv",nmarker = 1000,startn = 0)

}
