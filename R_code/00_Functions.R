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
cstats.table_PAM <- function(dist, k) {
  clust.assess <- c("cluster.number","wb.ratio","dunn","avg.silwidth")
  output.stats <- matrix(ncol = length(clust.assess), nrow = k-1)

   for(i in 2:k)
  {
    pam_fit <- pam(dist, diss = TRUE, k = i)
    pam.clust.num <- (pam_fit$clustering)

    output.stats[i-1,] <- unlist(cluster.stats(d = dist, clustering = pam.clust.num)[clust.assess])
  }
  
  colnames(output.stats) <- clust.assess
  rownames(output.stats) <- c(2:k)

  #return(output.stats)
  output.stats.df <- as.data.frame(output.stats)
  
  resultado <- median(c(output.stats.df[which.max(output.stats.df$avg.silwidth),]$cluster.number,
                        output.stats.df[which.max(output.stats.df$dunn),]$cluster.number,
                        output.stats.df[which.min(output.stats.df$wb.ratio),]$cluster.number))
  
  return(resultado)
}
