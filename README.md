# ClustALL_AD
Code repository for the analysis in the following article:

###### ClustALL: A robust clustering strategy for stratification of patients with acutely decompensated cirrhosis.

#### Medarxive link:


- [1. Context and Aim](##1.-context-and-aim.)  
- [2. Programming Environment & Scripts](##2.-programming-environment-&-scripts.)  
- [3. Input data](##3.-input-data.)  
- [4. Approach](##4.-approach.)  
- [5. Data used](##5.data-used.)  

## 1. Context and Aim.
Patient heterogeneity poses a significant obstacle in both individual patient management and the design of clinical trials, particularly in the management of complex diseases. Many existing clinical classifications rely on scores constructed for predicting patient outcomes. However, these conventional methods may overlook features contributing to heterogeneity that do not necessarily translate into prognostic implications.

To tackle the issue of patient heterogeneity upon hospital admission considering clinical features, we introduced ClustALL, a computational pipeline adept at handling common challenges associated with clinical data, such as mixed data types, missing values, and collinearity. ClustALL also facilitates the unsupervised identification of multiple and robust stratifications. 

 ![alt text](https://github.com/TranslationalBioinformaticsUnit/ClustALL_AD/blob/main/Fig_1_6_25.jpg?raw=true)



## 2. Programming Environment & Scripts.
ClustALL pipeline has been stablished in R version 4.2.2 (2022-10-31).

- 00_Functions.R -> Contains all the internal functions and required packages that are necessary to run the ClustALL pipeline.
- 00_Imputation.R -> Contains an example on performing imputation with MICE package
- 01a_ClustALL_noImputation.R -> Complete cluster pipeline for data sets with **NO missing values**. 
   >**Step 1**: *Data Complexity Reduction*
   >**Step 2**: *Stratification Process*
   >**Step 3**: *Stratification Representatives*.
   
- 01b_ClustALL_Imputation.R -> Complete cluster pipeline for data sets with **missing values**. The pipeline performs **imputations** to deal with NAs (the number of imputations can be specified). 
   >**Step 0:** *Imputation*
   >**Step 1**: *Data Complexity Reduction*
   >**Step 2**: *Stratification Process*
   >**Step 3**: *Stratification Representatives*.

## 3. Input data.
ClustALL is capable of processing both binary and numerical clinical variables as its input. Categorical features undergo transformation using a one-hot encoder method. A minimum of two features is necessary, although incorporating additional features enhances the precision of clustering. It is crucial to acknowledge that augmenting the number of features might also escalate computation time.

e.g.   

 ![alt text](https://github.com/TranslationalBioinformaticsUnit/ClustALL_AD/blob/main/Fig_2.jpg?raw=true)
 
## 4. Approach. 
The ClustALL methodology consists of three main steps: 

*(1) Data Complexity Reduction* (depicted in the Green Panel in the Figure) aims to simplify the original dataset by reducing redundant information, specifically highly correlated variables. This step yields a series of embeddings, each derived from distinct groupings of clinical variables.   

*(2) Stratification Process* (depicted in the Purple Panel of the Figure) involves multiple stratification analyses for each embedding. Various combinations of distance metrics and clustering methodologies are utilized in this step, generating stratifications denoted as "embedding + distance metric + clustering method".  

*(3) Consensus-based Stratifications* (depicted in the Red Panel of the Figure) aims to identify robust stratifications with minimal variation when parameters ("embedding + distance metric + clustering method") are slightly adjusted. ClustALL conducts a population-based robustness analysis for each stratification using bootstrapping. Non-robust combinations are excluded, and resulting stratifications are compared using the Jaccard distance. A heatmap visually highlights representative stratifications (indicated by green squared lines). The selection of representative stratifications ensures the preservation of those demonstrating parameter-based robustness even when various parameters, such as distance metrics or clustering methods, are altered. A final stratification for each group is chosen based on the centroid (green squares).

## 5. Data used.
As a proof of concept, CLustALL was applied to a subset of patients from the PREDICT (NCT03056612) cohort, a prospective European multicenter cohort comprising patients with acutely decompensated cirrhosis (AD). A total of 74 clinical variables from 766 patients collected at hospital inclusion with less than 30% of missing values were considered as input. 

Researchers who provide a methodology sound proposal can apply for the data, as far as the proposal is in line with the research consented by the patients. These proposals should be requested through www.datahub.clifresearch.com. Data requestors will need to sign a data transfer agreement. 


