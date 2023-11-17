### 01_ClustALL_noImputation.R
# This file contains the code to run ClustALL taking as an input a complete dataset with no missing values.

set.seed(1234)

######## Step 0. Prepare the environment ########
# Load the packages and functions described in the "00_Functions.R" file
source("00_Functions.R")

# Load the complete dataset
load("../data/data_use.RData") 

# Dataset exploration
print("Dataset dimension -") 
dim(data_use)

print("Position of missing values -") # Data frame should not contain missig values. If it has missing values, go to script 01b.
which(is.na(data_use)) 

print("Variables type -")
str(data_use)

nvariables <- ncol(data_use) # Store the number of variables

# Categorical variables are converted to binary with one-hot enconding
tobin <- data_use[, sapply(data_use, is.factor)] # Select factor variables
# Define one-hot encoding function
dummy <- dummyVars(" ~ .", data=tobin)
# Perform one-hot encoding on data frame
bin <- data.frame(predict(dummy, newdata=tobin))
names(bin) <- names(tobin)

data_use <- data_use[ , !(names(data_use) %in% names(tobin))] # Remove factor variables
data_use <- cbind(data_use, bin) # Add binarized facor variables

data_use[,1:nvariables] <- data.frame(apply(data_use[,1:nvariables], 2, as.numeric))

# Create the objects to store the ClustALL stratifications (a: Correlation + K-means ,b: Correlation + H-clust, c: Gower +  K-medoids, d: Gower + H-clust)
summary_clusters_a <- matrix(0,nvariables,1)
summary_clusters_b <- matrix(0,nvariables,1)
summary_clusters_c <- matrix(0,nvariables,1)
summary_clusters_d <- matrix(0,nvariables,1)
​
summary_matrices_a <- vector("list",nvariables)
summary_matrices_b <- vector("list",nvariables)
summary_matrices_c <- vector("list",nvariables)
summary_matrices_d <- vector("list",nvariables)
​
for(i in 1:(nvariables))
{
  m1 <- matrix(0,nrow(data_use),nrow(data_use))
  rownames(m1) <- rownames(data_use)
  colnames(m1) <- rownames(data_use)
  
  m2 <- matrix(0,nrow(data_use),nrow(data_use))
  rownames(m2) <- rownames(data_use)
  colnames(m2) <- rownames(data_use)
  
  m3 <- matrix(0,nrow(data_use),nrow(data_use))
  rownames(m3) <- rownames(data_use)
  colnames(m3) <- rownames(data_use)
  
  m4 <- matrix(0,nrow(data_use),nrow(data_use))
  rownames(m4) <- rownames(data_use)
  colnames(m4) <- rownames(data_use)
  
  summary_matrices_a[[i]] <- m1
  summary_matrices_b[[i]] <- m2
  summary_matrices_c[[i]] <- m3
  summary_matrices_d[[i]] <- m4
}
​
######## Step 1. Data Complexity Reduction ########
### Step 1.1. Dendrogram to assess the correlation between variables, and stablish the number of depths

cor_est_var <- cor(data_use, method="spearman") # Obtain correlations between variables
variables_clust  <- hclust(as.dist(cor_est_var)) # Compute the dendrogram for variables by hierarchical clustering
possible_heights <- variables_clust$height # Obtain the number of depths, based on how many times the dendrogram can be cutted

### Step 1.2. Reprocessing
for (heights_cut in (length(possible_heights)-(ncol(data_use)-2)):(length(possible_heights)-1)) { # compute for all possible depths
  treego <- cutree(variables_clust, h=possible_heights[heights_cut]) 
  groups_go <- unique(treego)
  variables_use <- names(treego[treego==groups_go[1]]) # choose first set of variables
  
  if (length(variables_use)==1) { # not generate PCA with one variable, if there is one variable alone, keep it as it is
    data_PCA <- data_use[,variables_use]
    
  } else { # generate the PCA 
    data_PCA <- generatePCA_derived(data=data_use[,variables_use],
                                    variability=50, maxvar=3)
  }
  
  if (length(groups_go) > 1) {
    for (h in 2:length(groups_go)) { # for every set of the tree
      variables_use <- names(treego[treego==groups_go[h]])
      
      if (length(variables_use)==1) { 
        data_PCA <- cbind(data_PCA, data_use[,variables_use])
        
      } else {
        data_PCA <- cbind(data_PCA, generatePCA_derived(data=data_use[,variables_use],
                                                        variability=50,maxvar=3))
      }
    } 
  }
  ## END OF Step 1. 
  
######## Step 2. Statrification Process ########
### Step 2.1. Conducting statrifications 
data_PCA_scaled <- apply(data_PCA, 2, scale) # scale PCA values for correlation distances; the clustering methods that are going to be applied here are (a) K-MEANS and (b) H-CLUST clustering
rownames(data_PCA_scaled) <- 1:nrow(data_PCA_scaled) # rename otherwise clValid function does not detect rownames of the df
  
#### Step 2.1.1. a) Correlation Distance +  K-Means clustering method
data_PCA.clValid_internal_kmeans <- clValid(data_PCA_scaled, 2:6, clMethods=c("kmeans"), 
                                              validation="internal",maxitems=nrow(data_PCA_scaled),metric="correlation") # we assess from 2 to 6 clusters
oS_kmeans <- as.numeric(mlv(optimalScores(data_PCA.clValid_internal_kmeans)[,3], method = "mfv"))
  
if (length(oS_kmeans) > 1) {
  oS_kmeans <- median(oS_kmeans) # consensus to select the best cluster number
}
  
kmeans_res <- clusters(data_PCA.clValid_internal_kmeans, "kmeans")
kmeans_res_c <- kmeans_res[as.character(oS_kmeans)][[1]]$cluster
  
for (t in 1:as.numeric(oS_kmeans)) { # for every k 
  induse <- as.numeric(names(kmeans_res_c[kmeans_res_c==t]))
  summary_matrices_a[[heights_cut]][induse, induse] <- summary_matrices_a[[heights_cut]][induse, induse] + 1 # add result
}
summary_clusters_a[heights_cut, 1] <- oS_kmeans
  
  
#### Step 2.1.2. b)  Correlation Distance + H-Clust clustering method
data_PCA.clValid_internal_hclust <- clValid(data_PCA_scaled, 2:6, clMethods=c("hierarchical"), 
                                              validation="internal",maxitems=nrow(data_PCA_scaled),metric="correlation")
oS_hclust <- as.numeric(mlv(optimalScores(data_PCA.clValid_internal_hclust)[,3], method = "mfv"))
  
if (length(oS_hclust) > 1) {
  oS_hclust <- median(oS_hclust) # consensus to select the best cluster number
} 
  
hclust_res <- clusters(data_PCA.clValid_internal_hclust, "hierarchical")
hclust_res_c <- cutree(hclust_res, k=oS_hclust)
  
for (t in 1:oS_hclust) {
  induse <- as.numeric(names(hclust_res_c[hclust_res_c==t]))
  summary_matrices_b[[heights_cut]][induse, induse] <- summary_matrices_b[[heights_cut]][induse, induse] + 1 # add result
}
summary_clusters_b[heights_cut, 1] <- oS_hclust
  
#### Gower Distance preparation
# Gower distances: the clustering methods that are going to be applied here are (c) K-Medoids and (d) H-Clust clustering
# Gower treats different the different data types, let's assess the type of the embeddings and perform the necessary modifications
PCA_gower <- as.data.frame(data_PCA) 
colnames(PCA_gower) <- c(1:dim(PCA_gower)[2])
binary <- which(apply(PCA_gower,2,function(x) {all(x %in% 0:1)})) # those variables with binary values
  
if (length(binary)==1) {
  PCA_gower[,binary] <- factor(PCA_gower[,binary])
    
} else {
  PCA_gower[,binary] <- data.frame(apply(PCA_gower[,binary], 2, as.logical))
}
  
PCA_gower <- PCA_gower %>% mutate_if(is.character, as.factor)
gower_dist <- daisy(PCA_gower, metric="gower")
  
  
#### Step 2.1.3. c) Gower Distance + K-Medoids clustering method
# apply the previous formula to get the optimal number of clusters for each depth of the dendrogram with gower distance + PAM (K-medoids)   
summary_clusters_c[heights_cut, 1] <- cstats.table_PAM(gower_dist, 6)
pam_res <- pam(gower_dist, k=summary_clusters_c[heights_cut, 1])
pam_res_c <- pam_res$clustering
  
# calculate for each depth and its corresponding optimal cluster number
for (t in 1:summary_clusters_c[heights_cut,  1]) {
  induse <- as.numeric(names(pam_res_c[pam_res_c==t]))
  summary_matrices_c[[heights_cut]][induse, induse] <- summary_matrices_c[[heights_cut]][induse, induse] + 1 # add result
}
  
  
#### Step 2.1.4. d) Gower Distance + H-Clust method 
divisive.clust <- diana(as.matrix(gower_dist), diss = TRUE, keep.diss = TRUE)
summary_clusters_d[heights_cut, 1] <- cstats.table_hclust(gower_dist, divisive.clust, 6)
hclustgow_res_c <- cutree(divisive.clust, k=summary_clusters_d[heights_cut, 1])
names(hclustgow_res_c) <- 1:dim(data_use)[1]
  
for (t in 1:summary_clusters_d[heights_cut, 1])  {
  induse <- as.numeric(names(hclustgow_res_c[hclustgow_res_c==t]))
  summary_matrices_d[[heights_cut]][induse, induse] <- summary_matrices_d[[heights_cut]][induse, induse] + 1 # add result
}
}

### END OF Step 2.1


# remove those rows with only 0s 
nclust_a <- as.data.frame(summary_clusters_a[which(apply(summary_clusters_a, 1, sum)>0),])
nclust_b <- as.data.frame(summary_clusters_b[which(apply(summary_clusters_b, 1, sum)>0),])
nclust_c <- as.data.frame(summary_clusters_c[which(apply(summary_clusters_c, 1, sum)>0),]) 
nclust_d <- as.data.frame(summary_clusters_d[which(apply(summary_clusters_d, 1, sum)>0),])

# put all the clusters together
summary_n_clust <- c(apply(nclust_a, 1, function(x) median(x, na.rm=T)),
                     apply(nclust_b, 1, function(x) median(x, na.rm=T)),
                     apply(nclust_c, 1, function(x) median(x, na.rm=T)),
                     apply(nclust_d, 1, function(x) median(x, na.rm=T)))

# renaming the matrices
names(summary_matrices_a) <- paste("cuts_a_",1:length(summary_matrices_a), sep="")
names(summary_matrices_b) <- paste("cuts_b_",1:length(summary_matrices_b), sep="")
names(summary_matrices_c) <- paste("cuts_c_",1:length(summary_matrices_c), sep="")
names(summary_matrices_d) <- paste("cuts_d_",1:length(summary_matrices_d), sep="")

summary_matrices_MEASURES <- c(summary_matrices_a, summary_matrices_b,
                               summary_matrices_c, summary_matrices_d)

# remove those matrices with only 0s
tmp <- lapply(summary_matrices_MEASURES, function(x) as.matrix(x[]))
if (length(which(lapply(tmp, sum)==0)) > 0) {
  summary_matrices_MEASURES = summary_matrices_MEASURES[-which(lapply(summary_matrices_MEASURES, function(x) sum(x[]))==0)]
  message("Removing the following empty matrices: ") 
  message(paste(names((which(lapply(tmp, sum)==0)))), ".")
}


######## Step 2.2. Filtering non-robust stratifications ########
summary_clusters <- vector("list", length(summary_matrices_MEASURES))

for(i in 1:length(summary_clusters)) {
  hclustgo <- hclust(1-as.dist(summary_matrices_MEASURES[[i]][]))
  summary_clusters[[i]] <- cutree(hclustgo, k=summary_n_clust[i])
}
names(summary_clusters) <- names(summary_matrices_MEASURES)

# Jaccard Distance
JACCARD_DISTANCE <- matrix(NA, length(summary_matrices_MEASURES), length(summary_matrices_MEASURES))
rownames(JACCARD_DISTANCE) <- names(summary_matrices_MEASURES)
colnames(JACCARD_DISTANCE) <- rownames(JACCARD_DISTANCE)

# For two clusterings of the same data set, this function calculates the similarity statistic specified of the clusterings from the comemberships of the observations. 
# Basically, the comembership is defined as the pairs of observations that are clustered together.
for(i in 1:nrow(JACCARD_DISTANCE)) {
  for(j in 1:nrow(JACCARD_DISTANCE)) {
    JACCARD_DISTANCE[i,j] <- cluster_similarity(summary_clusters[[i]], summary_clusters[[j]],
                                                similarity="jaccard", method ="independence")
  }
}

# Robustness of the clustering
summary_matrices_STABILITY <- matrix(NA, length(summary_matrices_MEASURES), 3)
for(i in 1:length(summary_matrices_MEASURES)){
  invisible(capture.output( # avoid printing function messages
    r1 <- clusterboot(data=as.dist(1000-summary_matrices_MEASURES[[i]]),
                      B=100, distances=TRUE, bootmethod="boot",
                      bscompare=TRUE, multipleboot=FALSE,
                      jittertuning=0.05, noisetuning=c(0.05,4),
                      subtuning=floor(nrow(data)/2),
                      clustermethod=disthclustCBI,noisemethod=FALSE,count=TRUE,
                      showplots=FALSE,dissolution=0.5,
                      recover=0.75,method="complete",k=2)
  ))
  
  summary_matrices_STABILITY[i,1:2] <- as.numeric(r1$bootmean[1:2])
  summary_matrices_STABILITY[i,3] <- mean(summary_matrices_STABILITY[i,1:2])
}


quantileuse <- 0.85 # Stratifications with less than 85% stability are excluded
qgo <- quantile(summary_matrices_STABILITY[,3], quantileuse)
JACCARD_DISTANCE_F <- JACCARD_DISTANCE[summary_matrices_STABILITY[, 3] >= qgo,summary_matrices_STABILITY[, 3]>=qgo]
### END OF Step 2.2. 

######## Step 3. Stratification Representatives ########
maxgo <- max(JACCARD_DISTANCE_F[JACCARD_DISTANCE_F < 1])
ordergo <- hclust(1 - as.dist(JACCARD_DISTANCE_F))$order

# This plot contains the robust stratifications
corrplot(JACCARD_DISTANCE_F[ordergo, ordergo],is.corr=F,method="shade",
         cl.lim=c(min(JACCARD_DISTANCE),1),main="",order="original",
         tl.col="black", tl.srt=45,cl.pos="b",cl.length=4)

# obtain the definitive clusters
definitive_clusters = obtainDefCluster(JACCARD_DISTANCE_F[ordergo, ordergo])

# chose the representative cluster with minimum population in each cluster. As default 0.05 (5%)
chooseClusters(definitive_clusters, summary_clusters, 0.05)

## END OF Step 3. 



                                                                      
