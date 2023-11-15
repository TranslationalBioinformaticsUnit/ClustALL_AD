######## Load required packages and functions ########


### Required packages

x <- c("dplyr", "ggplot2", "cluster", "mclust", "clValid", "modeest", "mice",
       "FactoMineR", "MASS", "corrplot", "clusteval", "fpc", "haven", "naniar",
      "reshape2", "circlize", "labelled",  "Rtsne")

lapply(x, require, character.only = TRUE)

###  Required Functions

# > PCA derived information

generatePCA_derived <- function(data,variability,maxvar)
{
  pca_go = PCA(data, graph = FALSE,ncp=min(c(maxvar,ncol(data))))
  variab <- pca_go$eig[,3]
  whichgo <- min(which(variab>50))
  varigo <- min(c(whichgo,maxvar))
  res <- pca_go$ind$coord[,c(1:varigo)]
  res
}

# > 
