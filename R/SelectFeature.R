################################################################
# SelectFeature: feature selection before applying PREDE model
################################################################

SelectFeature <- function(mat,method = "cv",nmarker=1000, startn=0){

  if(method == "topvar"){ ## sort sd
    Var = rowVars(mat)
    difgene = order(Var,decreasing = T)[startn + 1:nmarker]
    index.select = rownames(mat)[difgene]
  }
  if(method == "random"){ ## random
    index.select = sample(rownames(mat),size = nmarker)
  }
  if(method == "cv"){ ## sort sd/mean
    mm = rowMeans(mat)
    vv = rowVars(mat)
    cv = sqrt(vv) / (mm + 1)
    cv[is.na(cv)] = 0
    ix = sort(cv, dec=TRUE, index=TRUE)$ix
    index.select = rownames(mat)[ix[startn + 1:nmarker]]
  }
  return(index.select)
}




