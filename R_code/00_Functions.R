### 00_Functions.R
# This file contains the required packages and functions to run ClustALL


######## Check if the required packages are installed ########
packages <- c("dplyr", "ggplot2", "cluster", "mclust", "clValid", "modeest", "mice", 
              "FactoMineR", "MASS", "corrplot", "clusteval", "fpc", "haven", "naniar",
              "reshape2", "circlize", "labelled",  "Rtsne")

new.packages = packages[!(packages %in% installed.packages()[,"Package"])]

install.packages(new.packages)

######## Load required packages and functions ########
### Required packages

packages <- c("dplyr", "ggplot2", "cluster", "mclust", "clValid", "modeest", "mice",
       "FactoMineR", "MASS", "corrplot", "clusteval", "fpc", "haven", "naniar",
      "reshape2", "circlize", "labelled",  "Rtsne")

lapply(packages, require, character.only = TRUE)

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

       # > Function to choose the best cluster number considering the Gower distance and the K-medoids clustering algorithm
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
  
  res <- median(c(output.stats.df[which.max(output.stats.df$avg.silwidth),]$cluster.number,
                        output.stats.df[which.max(output.stats.df$dunn),]$cluster.number,
                        output.stats.df[which.min(output.stats.df$wb.ratio),]$cluster.number))
  
  return(res)
}

       # > Function to choose the best cluster number considering the Gower distance and the H-clust clustering algorithm
cstats.table_hclust <- function(dist, tree,k) {
  clust.assess <- c("cluster.number","wb.ratio","dunn","avg.silwidth")
  output.stats <- matrix(ncol = length(clust.assess), nrow = k-1)

  for(i in 2:k)
  {
    
    output.stats[i-1,] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.assess])
   
      
  }
  
  colnames(output.stats) <- clust.assess
  rownames(output.stats) <- c(2:k)

  #return(output.stats)
  output.stats.df <- as.data.frame(output.stats)
  
  res <- median(c(output.stats.df[which.max(output.stats.df$avg.silwidth),]$cluster.number,
                        output.stats.df[which.max(output.stats.df$dunn),]$cluster.number,
                        output.stats.df[which.min(output.stats.df$wb.ratio),]$cluster.number))
  
  return(res)
}
