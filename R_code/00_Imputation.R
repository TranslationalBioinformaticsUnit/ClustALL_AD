### 01_ClustALL_noImputation.R
# This file prepares the data with missing values for imputation with MICE R package
# You should run ClustALL taking as an input the output of the imputation and consider the framework the 01b_ClustALL_Imputation to compute clustall.

# Note that imputation would depend specifically on the nature of your data. This script shows the imputation perform in the data mentioned in 05.Dataset
# ClustALL can take impute data perform in other platforms.

#########  Dealing with missing data: Multivariate imputation through chained equations (MICE) ######## 

### We first explored and imputed the data due to the presence of missing values. The possible dependence among variables was considered for the imputation.
## Raw data loading with missing data
data_to_impute <- readRDS("~/data_to_impute.rds")

# Analysis of the missing values:

table(is.na(data_to_impute)) #number of NAs in the dataset
colSums(is.na(data_to_impute)) # columns with NAs

# Remove columns with more than 30% NA
data_to_impute <- data_to_impute[, which(colMeans(!is.na(data_to_impute)) > 0.3)]

# Display missing-data patterns
md.pattern(data_to_impute)

# Modified the predictor matrix to select the predictor for each variable:

# Prepare the predictor matrix for the imputation:
# maximum number of iterations maxit set to zero, a fast way to create a mids object containing the default settings
ini <- mice(data_to_impute[,c(1:74)], maxit = 0, print = FALSE)

# Selection of predictors for each variable, if we know that one variable depends on other, we take it into account
pred <- ini$pred # predictor matrix

pred[c(1:45,47:51,53), ] <- 0 # remove variables without missing values of the matrix

pred["CHILD", -c(20:21,48:49)] <- 0
pred["SPO2", -c(38,59)] <- 0
pred["SPFIO2", -c(38,39,72)] <- 0
pred["bmi", -c(1,2,53,54)] <- 0
pred["LYM", -47] <- 0
pred["MONO", -47] <- 0
pred["NEUT", -47] <- 0
pred["PLAT", -c(47:51,59:63,65:71)] <- 0
pred["ALB", -c(29,30,17,43,48,60,66)] <- 0
pred["GGT", -c(18,25,27,37,57,58)] <- 0
pred["K", -c(1,2,17,44,48:51,64)] <- 0
pred["CRP", -c(1,2,24,56,70)] <- 0
pred[c(52,54,56:58,70,66), c(8:11,16,28,40,41,53:56)] <- 0 #ADHIST, tobacco, alcohol, active_alcohol, wgt, alt
pred["AST", -66] <- 0
pred["HGB", -60] <- 0
pred["HCT", -59] <- 0

## Criteria for performing ths specific imputation due to the nature of the variables:
  # CHILD depends on BILI, ALB, INR, HE and ASCI
  # SPO2 depends on hgb and RESPDYS
  # SPFIO2 depends on SPO2, RESPDYS and RESPFLR
  # BMI depends on hgt, wgt, age and sex
  # LYM depends on WBC
  # MONO depends on WBC
  # NEU depends on WBC
  # PLAT laboratory previous variables as predictors
  # ALB depends on LIVERDYS and LIVERFLR - ALB correlates with etiology6, MELDNA, INR, HCT and ALT
  # AST and ALT are related: predict ALT (less NAs) with the other variables and predict AST with ALT
  # GGT depends on alcohol, HCC, CARDIFLR, DM and ASH
  # K depends on Na - K correlates (very low) with age, sex, etiology6, CLIFAD, INR, BILI, CREAT and PLAT
  # CRP depends on age, sex, GLUC, tobacco and BINF

# imputation with 5 iterations and 1000 imputations
nimp <- 1000 # number of imputations
imp <- mice(data_to_impute, m=nimp, pred=pred, maxit=5, seed=1234, print=FALSE)

# save result
saveRDS(imp, file="imp.Rda")

# Check predictors for each variable:
imp$formulas

## Checking imputations
# Plot the densities for all variable with missing values
densityplot(imp, scales = list(x = list(relation = "free")), layout = c(5, 1))

# Plor for one variable
densityplot(imp, ~K, scales = list(x = list(relation = "free")))

# Stripplot for all variables in the dataset
stripplot(imp, col = c("grey", mdc(2)), pch = c(1, 20))

# Scatterplot (separeted observed in blue and imputed in red data)
xyplot(imp_pred, alcohol ~ sex | .imp, pch = c(1, 20), cex = c(0.8, 1.2), scales = list(tick.number = 3))
