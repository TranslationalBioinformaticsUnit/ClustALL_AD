# ClustALL_AD
ClustALL: A robust clustering strategy for stratification of patients with acutely decompensated cirrhosis.

## 1. Context and Aim.
Patient heterogeneity poses a significant obstacle in both individual patient management and the design of clinical trials, particularly in the management of complex diseases. Many existing clinical classifications rely on scores constructed for predicting patient outcomes. However, these conventional methods may overlook features contributing to heterogeneity that do not necessarily translate into prognostic implications.

To tackle the issue of patient heterogeneity upon hospital admission considering clinical features, we introduced ClustALL, a computational pipeline adept at handling common challenges associated with clinical data, such as mixed data types, missing values, and collinearity. ClustALL also facilitates the unsupervised identification of multiple and robust stratifications. 

 Image
 
### Medarxive link:

## 2. Programming Environment & Scripts.
ClustAll pipeline has been stablished in R version 4.2.2 (2022-10-31).

- 00_Functions.R -> Contains all the internal functions and required packages that are necessary to run the ClustAll pipeline.
- 01a_ClustAll_noImputation.R -> Complete cluster pipeline for data sets with **NO missing values**. 
   >**Step 1**: *Data Complexity Reduction*  **Step 2**: *Stratification Process*  **Step 3**: *Stratification Representatives*.
   
- 01b_ClustAll_Imputation.R -> Complete cluster pipeline for data sets with **missing values**. The pipeline performs **imputations** to deal with NAs (the number of imputations can be specified). 
   >**Step 0:** *Imputation*  **Step 1**: *Data Complexity Reduction*  **Step 2**: *Stratification Process*  **Step 3**: *Stratification Representatives*.

## 3. Input data.
ClustALL is capable of processing both binary and numerical clinical variables as its input. Categorical features undergo transformation using a one-hot encoder method. A minimum of two features is necessary, although incorporating additional features enhances the precision of clustering. It's crucial to acknowledge that augmenting the number of features might also escalate computation time.

## 4. Approach.
- 00_Functions.R -> Contains all the internal functions and required packages that are necessary to run the ClustAll pipeline.
- 01a_ClustAll_noImputation.R -> Complete cluster pipeline for data sets with **NO missing values**. 
   >**Step 1**: *Data Complexity Reduction*  **Step 2**: *Stratification Process*  **Step 3**: *Stratification Representatives*.
   
- 01b_ClustAll_Imputation.R -> Complete cluster pipeline for data sets with **missing values**. The pipeline performs **imputations** to deal with NAs (the number of imputations can be specified). 
   >**Step 0:** *Imputation*  **Step 1**: *Data Complexity Reduction*  **Step 2**: *Stratification Process*  **Step 3**: *Stratification Representatives*.
   
## 5. Data used.
As a proof of concept, CLustAll was applied to PREDICT [ref] a prospective European multicenter cohort comprising patients with acutely decompensated cirrhosis (AD) (n=766). A total of 74 clinical variables collected at hospital inclusion with less than 30% of missing values were considered as input. 

ClustAll revealed five robust stratifications based solely on clinical data collected at hospital admission. All identified stratifications incorporated markers of impaired liver function and the number of organ dysfunction or failure, with most also encompassing precipitating events. Upon focusing on one of these stratifications, patients were categorized into three clusters characterized by typical clinical features, which also demonstrated prognostic value. A re-assessment of patient stratification during follow-up delineated outcomes, further enhancing the prognostic value of the stratification. These findings were also validated in ACLARA [ref], an independent prospective multicenter cohort of patients from Latin America (n=580).



