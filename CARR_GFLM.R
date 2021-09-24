CARR_GFLM <- function(X,      # data.frame in long format with five labeled columns (described below)
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
                      a,      # covariate domain (vector)
                      Basisw, # functional basis (matrix)
                      Pw,     # functional penalty (matrix)
                      Basisa, # covariate basis (matrix)
                      Pa,     # covariate penalty (matrix)
                      Basisr, # regional basis (matrix)
                      Pr,     # regional penalty (matrix)
                      dist,   # distribution of outcome (character)
                      pvar,   # proportion of variation (scalar in [0, 1])
                      ci,     # return confidence intervals (1 = yes, 0 = no)
                      kbasis, # degrees of freedom for region-referenced mean surface smoothing, c(functional, covariate)
                      weights # implement weighted smoothing for region-referenced mean surface (1 = yes, 0 = no), for binomial outcomes only
                      ){ 
#############################################################################
## Description: fits covariate-adjusted region-referenced generalized functional linear model (CARR-GFLM) 
## described in "Covariate-Adjusted Region-Referenced Generalized Functional Linear Model for EEG Data"
## by Scheffler et al. (2019).
## Inputs:      see above
## Returns:     list() 
##              rm: mgcv model fit (list)
##              VCov: covariance matrix for model parameters (matrix)
##              V: right singular vectors of design matrix (matrix)
##              fit: fitted values returned from mgcv model (vector)
##              betaHat: regression function (list)
##              betaUL: upper bound for Bayesian 95% point-wise confidence interval (list)
##              betaLL: lower bound for Bayesian 95% point-wise confidence interval (list)
##              regMod: region-referenced smooth models (list)
##              regMean: region-referenced smooths (list)
## CARR-GFLM  Outline:
##              0. Calculate summary variables and format data
##              1. Mean center subject-specific functional predictors
##              2. Construct design matrix
##              3. Construct penalty matrix
##              4. Fit CARR-GFLM using mgcv
##              5. Form regression function with Bayesian 95% point-wise confidence intervals
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
    wts[y ==as.numeric(names(table(y))[as.numeric(which.min(table(y)))])] <- max(1/(table(y)[1]/table(y)[2]), 1/(table(y)[2]/table(y)[1]))
  }
}
wts <- rep(wts, each = lw)

#############################################################################
# 1. Mean center subject-specific functional predictors
#############################################################################
print("1. Mean center subject-specific functional predictors...")

regMean <- matrix(list(), nrow = R)
regMod <- matrix(list(), nrow = R)
newData <- data.frame(cbind(rep(w, length(a)), rep(a, each = length(w))))
colnames(newData) <- c("func", "covar")
for(r in 1:R){
  regMod[[r]] <- gam(x ~ te(func, covar, k = kbasis, bs = "ps", m = 2), data = X[which(X$reg == r), ], weights = wts)
  X$x[which(X$reg == r)] <- X$x[which(X$reg == r)] - regMod[[r]]$fitted.values
  regMean[[r]] <- matrix(predict(regMod[[r]], newdata = newData), nrow = lw)
}

Xmat <- t(matrix(X$x, nrow = lw * R)) # form matrix with subject-specific functional predictor in each row

#############################################################################
# 2. Construct design matrix
#############################################################################
print("2. Construct design matrix...")

# Form design matrix
B <- kron(Basisr, Basisw) # form partial tensor basis from regional and functional basis functions 
L <- matrix(NA, nrow = n, ncol = dim(Basisa)[2]) # extract evaluations of covariate basis functions at a_i 
index <- match(unlist(ai[, 2]), a)              
for(i in 1:n){
  L[i, ] <- Basisa[index[i], ]
}
Z <- by * Xmat %*% B # numerical integration of functional predictor and partial tensor basis
Zc <- matrix(NA, nrow = dim(Xmat)[1], ncol = ncol(Z) * ncol(Basisa))
for(i in 1:n){       # row kronecker product with covariate basis functions at a_i
  Zc[i, ] <- kron(L[i, ], Z[i, ])
}

# Reduce dimension of design matrix via SVD
svdZc <- svd(Zc) # perform SVD
V <- svdZc$v[, cumsum(svdZc$d^2 / sum(svdZc$d^2)) < pvar] # extract right singular vectors 
ZVc <- Zc %*% V  # reduced dimension design matrix
modDat <- list(ZVc = ZVc) # store design matrix in list 


#############################################################################
# 3. Construct penalty matrix
#############################################################################
print("3. Construct penalty matrix...")

# Penalty matrices in kronecker sum
P1 <- kron(Pa, kron(eye(dim(Basisr)[2]), eye(dim(Basisw)[2]))) 
P2 <- kron(eye(dim(Basisa)[2]), kron(Pr, eye(dim(Basisw)[2])))
P3 <- kron(eye(dim(Basisa)[2]), kron(eye(dim(Basisr)[2]), Pw)) 

# Adjust penalty matrices to account for dimension reduction via SVD
P1V <- t(V) %*% P1 %*% V 
P2V <- t(V) %*% P2 %*% V 
P3V <- t(V) %*% P3 %*% V

PP <- list(ZVc = list(P1V, P2V, P3V)) # store penalty matrix in list

#############################################################################
# 4. Fit CARR-GFLM using mgcv
#############################################################################
print("4. Fit CARR-GFLM using mgcv...")

rm <- gam(y ~ 1 + ZVc, data = modDat, paraPen = PP, family = dist, method = "REML")
VCov <- vcov(rm, freq = FALSE)[-1, -1] # store covariance matrix for model parameters (excluding intercept)

#############################################################################
# 5. Form regression function with Bayesian point-wise confidence intervals
#############################################################################
print("5. Construct the regression function with Bayesian point-wise confidence intervals...")

coefs <- c(unlist(rm$coefficients))[-1] # extract model coefficients (excluding intercept)
Tbasis = kron(Basisa, kron(Basisr, Basisw)) # construct tensor basis
betaTemp <- matrix((Tbasis %*% V %*% (coefs)), nrow = lw) # form vectorized regression function
betaHat <-  matrix(list(), nrow = R) # organize regression function by regions
for(r in 1:R){
  betaHat[[r]] <- matrix(0, nrow = lw, ncol = length(a))
  for(j in 1:length(a)){
    betaHat[[r]][, j] <- betaTemp[ ,((j - 1) * R + r)] 
  }
}


if(ci == 1){ # calculate Bayesian point-wise confidence intervals
  betaVar <- rep(NA, dim(Tbasis)[1]) # calculate point-wise variance of regression function
  for(j in 1:length(betaVar)){ # vectorized Bayesian point-wise variance
    betaVar[j] <- Tbasis[j, ] %*% V %*% VCov %*% t(V) %*% Tbasis[j, ] 
  }
  temp <- matrix(betaVar, nrow = length(w))
  betaVarMat <- matrix(list(), nrow = R)
  betaUL <- matrix(list(), nrow = R)
  betaLL <- matrix(list(), nrow = R)
  for(r in 1:R){ # organize point-wise variance by regions
    betaVarMat[[r]] <- matrix(0, nrow = length(w), ncol = length(a))
    for(j in 1:length(a)){
      betaVarMat[[r]][, j] <- temp[ ,((j - 1) * R + r)]
    }
    betaUL[[r]] <- betaHat[[r]] + 1.96 * sqrt(betaVarMat[[r]]) # upper 95% interval
    betaLL[[r]] <- betaHat[[r]] - 1.96 * sqrt(betaVarMat[[r]]) # lower 95% interval
  }
} else{
  betaUL <- NA
  betaLL <- NA
}
  
  
output <- list(rm,                # CARR-GFLM model (list)
               VCov,              # covariance matrix for tensor basis coefficients (matrix)
               V,                 # right singular vectors explaining 95% of the total variation (matrix)
               rm$fitted.values,  # fitted values (vector)
               betaHat,           # regression function (list)
               betaUL,            # regression function Bayesian confidence interval upper limit (list)
               betaLL,            # regression function Bayesian confidence interval lower limit (list)
               regMod,            # region-referenced smooth models (list)
               regMean            # region-referenced smooths (list)
               ) 
names(output) <- c("rm", "VCov", "V", "fit", "betaHat", "betaUL", "betaLL", "regMod", "regMean")
return(output)

}
