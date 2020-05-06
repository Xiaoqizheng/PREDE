pkgname <- "PREDE"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('PREDE')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("GenerateMixSample")
### * GenerateMixSample

flush(stderr()); flush(stdout())

### Name: generate_bulk
### Title: generate the mixture samples
### Aliases: generate_bulk

### ** Examples


## load example data
data(lung_exp)
W <- lung_exp[,1:6]

## generate the mixed samples based on the profile matrix of the cell types W
bulk <- generate_bulk(W,nSample =100,csd = 0.1)




cleanEx()
nameEx("GetCelltypeNum")
### * GetCelltypeNum

flush(stderr()); flush(stdout())

### Name: GetCelltypeNum 
### Title: get optimal number of total cell types in the mixture tumor
###   sample
### Aliases: 'GetCelltypeNum '

### ** Examples


## load example data
data(lung_exp)
W <- lung_exp[,1:6]

## generate the mixed samples
bulk <- generate_bulk(W,nSample =100,csd = 0.1)

## select the feature
feat <- select_feature(mat = bulk$Y,method = "cv",nmarker = 1000,startn = 0)

## get optimal number of total cell types
OptimalK <- GetCelltypeNum(bulk$Y[feat,],W=NULL,W1=W[feat,1:4],maxK = 10)



cleanEx()
nameEx("PREDE")
### * PREDE

flush(stderr()); flush(stdout())

### Name: PREDE
### Title: Partial-reference based deconvolution model
### Aliases: PREDE

### ** Examples


## load example data
data(lung_exp)
W <- lung_exp[,1:6]

## Partial-reference based deconvolution without total cell types W
bulk <- generate_bulk(W,nSample =100,csd = 0.1)
feat <- select_feature(mat = bulk$Y,method = "cv",nmarker = 1000,startn = 0)
PREDE(bulk$Y[feat,],W1=W[feat,1:4],type = "GE",K=7,iters = 100,rssDiffStop=1e-5)




cleanEx()
nameEx("SelectFeature")
### * SelectFeature

flush(stderr()); flush(stdout())

### Name: select_feature
### Title: Select feature for profile matrix of mixed samples
### Aliases: select_feature

### ** Examples


## load example data
data(lung_exp)
W <- lung_exp[,1:6]

## Partial-reference based deconvolution without total cell types W
bulk <- generate_bulk(W,nSample =100,csd = 0.1)
feat = select_feature(mat = bulk$Y,method = "cv",nmarker = 1000,startn = 0)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
