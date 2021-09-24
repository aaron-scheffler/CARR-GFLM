mGFLMi <- function(X,        # data.frame in long format with five labeled columns (described below)
                             # and row length equal to the length of the vectorized region-referenced
                             # functional predictors across all subjects
                                          # DATA.FRAME COLUMNS:
                                          # ID: subject ID 
                                          # x: region-referenced functional predictor
                                          # reg: regional argument
                                          # func: functional argument
                                          # covar: covariate argument
                      y,      # outcome (vector)
                      w,      # functional domain (vector)
                      Basisw, # functional basis (matrix)
                      Pw,     # functional penalty (matrix)
                      dist,   # distribution of outcome (character)
                      pvar,   # proportion of variation (scalar in [0, 1])
                      kbasis, # degrees of freedom for region-referenced mean surface smoothing, c(functional)
                      weights # implement weighted smoothing for region-referenced mean surface (1 = yes, 0 = no), for binomial outcomes only
                      ){ 
#############################################################################
## Description: fits multivariate generalized functional linear model with age interaction (m-GFLMi)
## described in "Covariate-Adjusted Region-Referenced Generalized Functional Linear Model for EEG Data"
## by Scheffler et al. (2019).
## Inputs:      see above
## Returns:     list() 
##              rm: mgcv model fit (list)
##              VCov: covariance matrix for model parameters (matrix)
##              V: right singular vectors of design matrix (matrix)
##              fit: fitted values returned from mgcv model (vector)
##              betaHat: regression function (matrix)
##              betaHati: regression function for interaction term (matrix)  
##              regMod: region-referenced smooth models (list)
##              regMean: region-referenced smooths (list)
## m-GFLMi  Outline:
##              0. Calculate summary variables and format data
##              1. Mean center subject-specific functional predictors
##              2. Construct design matrix
##              3. Construct penalty matrix
##              4. Fit m-GFLMi using mgcv
##              5. Construct regression functions
#############################################################################

# Install missing packages
list.of.packages <- c("mgcv", "splines", "data.table", "ppls", "pracma")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) 
    
# Load packages  
library(data.table)
library(splines)
library(mgcv)
library(ppls)
library(pracma)
library(Matrix)
  
#############################################################################
# 0. Calculate summary variables and format data
#############################################################################
print("O. Calculate summary variables and format data...")

# Summary variables
lw <- length(w) # length of functional domain
by <- w[2] - w[1] # width of funcitonal grid points
n <- length(unique(X$ID)) # sample size
R <- length(unique(X$reg)) # total regions 
ai <- unique(X[, c("ID", "covar")]) # covariate values by subject
ai[, 2] <- ai[order(ai[, 1]), 2]  # order covariate values by subject

# Data formatting
y = y[order(unique(X$ID))]
X = data.table(X)  # convert data matrix to data.table format 
X = X[order(ID, reg, func)] # order by subject, region, frequency

# Weights for smoothing by group
wts <- rep(1, n)
if(dist == "binomial"){
  if(weights == 1){
    wts[y == as.numeric(names(table(y))[as.numeric(which.min(table(y)))])] <- max(1/(table(y)[1]/table(y)[2]), 1/(table(y)[2]/table(y)[1]))
  }
}

wts <- rep(wts, each = lw)
#############################################################################
# 1. Mean center subject-specific functional predictors
#############################################################################
print("1. Mean center subject-specific functional predictors...")

regMean <- matrix(list(), nrow = R)
regMod <- matrix(list(), nrow = R)
newData <- data.frame(w)
colnames(newData) <- c("func")
for(r in 1:R){
  regMod[[r]] <- gam(x ~ te(func, k = kbasis, bs = "ps", m = 2), data = X[which(X$reg == r), ], weights = wts)
  X$x[which(X$reg == r)] <- X$x[which(X$reg == r)] - regMod[[r]]$fitted.values
  regMean[[r]] <- matrix(predict(regMod[[r]], newdata = newData), nrow = lw)
}

#############################################################################
# 2. Construct design matrix
#############################################################################
print("2. Construct design matrix...")

# Form basis design matrix
B <- Basisw # form basis

# Form design matrix for each region
Xmat <- matrix(list(), nrow = R) # subject-specific functional predictors
Zc <- matrix(list(), nrow = R)   # functional main effects design matrix
Zci <- matrix(list(), nrow = R)  # functional main effects*covariate design matrix

for(r in 1:R){
Xmat[[r]] <- t(matrix(X$x[which(X$reg == r)], nrow = lw)) # form matrix with subject-specific functional predictor in each row
Zc[[r]] <- by * Xmat[[r]] %*% B # numerical integration of functional predictor and basis
Zci[[r]] <- by * Xmat[[r]] %*% B # numerical integration of functional predictor and basis (for interaction term)
  for(i in 1:n){       # form interaction term between subject-specific functional predictor and age
    Zci[[r]][i, ] <- Zci[[r]][i, ] * as.numeric(ai[i,2])
  }

}


# Form an overall design matrix
Zc <- matrix(unlist(Zc), nrow = n)
Zci <- matrix(unlist(Zci), nrow = n)
Zc <- cbind(Zc, Zci) # concatenate main effects and interaction term


# Reduce dimension of design matrix via SVD
svdZc <- svd(Zc) # perform SVD
V <- svdZc$v[, cumsum(svdZc$d^2 / sum(svdZc$d^2)) < pvar] # extract right singular vectors 
ZVc <- Zc %*% V  # reduced dimension design matrix
modDat <- list(ZVc = ZVc, a = as.numeric(unlist(ai[,2]))) # store design matrix in list 

#############################################################################
# 3. Construct penalty matrix
#############################################################################
print("3. Construct penalty matrix...")

PP <- matrix(list(), nrow = R)  # main effect penalty matrix
PPi <- matrix(list(), nrow = R) # interaction term penalty matrix
PV <- matrix(list(), nrow = R)  # main effect penalty matrix adjusted for SVD
PVi <- matrix(list(), nrow = R) # interaction term penalty matrix adjusted for SVD
db <- dim(Basisw)[2]

for(r in 1:R){
  PP[[r]] <- matrix(0, nrow = db * R, ncol = db * R) # form penalty matrix
  PP[[r]][(1 + db * (r - 1)):((db) + db * (r - 1)), (1 + db * (r - 1)):((db) + db * (r - 1))] <- Pw 
  PP[[r]] <- bdiag(PP[[r]], matrix(0, nrow = db * R, ncol = db * R))
  PV[[r]] <- as.matrix(t(V) %*% PP[[r]] %*% V) # adjust penalty matrices to account for dimension reduction via SVD
  
  PPi[[r]] <- matrix(0, nrow = db * R, ncol = db * R) # form penalty matrix
  PPi[[r]][(1 + db * (r - 1)):((db) + db * (r - 1)), (1 + db * (r - 1)):((db) + db * (r - 1))] <- Pw 
  PPi[[r]] <- bdiag(matrix(0, nrow = db * R, ncol = db * R), PPi[[r]])
  PVi[[r]] <- as.matrix(t(V) %*% PPi[[r]] %*% V) # adjust penalty matrices to account for dimension reduction via SVD
}

PP <- list(ZVc = c(PV, PVi)) # store penalty matrix in list

#############################################################################
# 4. Fit CARR-GFLM using mgcv
#############################################################################
print("4. Fit m-GFLMi using mgcv...")

rm <- gam(y ~ 1 + a + ZVc, data = modDat, paraPen = PP, family = dist, method = "REML")

#############################################################################
# 5. Construct the regression function
#############################################################################
print("5. Construct the regression function...")

coefs <- c(unlist(rm$coefficients))[-1] # extract model coefficients (excluding intercept)
Tbasis = Basisw # construct tensor basis
db <- dim(Basisw)[2]

betaHat <-  matrix(NA, nrow = length(w), ncol = R) # organize regression function by regions
for(r in 1:R){
  betaHat[, r] <- Tbasis %*% (V %*% coefs[-1])[(1 + db * (r - 1)):((db) + db * (r - 1))]
}


betaHati <-  matrix(NA, nrow = length(w), ncol = R) # organize regression function by regions
for(r in 1:R){
  betaHati[, r] <- Tbasis %*% (V %*% coefs[-1])[((1 + db * R) + db * (r - 1)):((db * R  + db) + db * (r - 1))]
}

  
output <- list(rm,               # m-GFLM model (list)
               rm$fitted.values, # fitted values (vector)
               betaHat,          # main effect regression function (matrix)
               betaHati,         # interaction term regression function (matrix)
               V,                # right singular vectors explaining 95% of the total variation (matrix)
               regMod,           # region-referenced smooth models (list)
               regMean           # region-referenced smooths (list)
               )


names(output) <- c("rm","fit", "betaHat", "betaHati", "V", "regMod", "regMean")
return(output)


}
