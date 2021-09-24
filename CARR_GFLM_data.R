CARR_GFLM_data <- function(n,       # sample size
                           rho,     # level of regional dependency 
                           R,       # number of regions
                           w,       # functional domain
                           a        # covariate domain
                           ){
  
#############################################################################
## Description: Function for simulating data for the GARR_GFLM.R function as described
##              in the simulation section found in Section 4 of "Covariate-Adjusted 
##              Region-Referenced Generalized Functional Linear Model for EEG Data" 
##              by Scheffler et al. (2018).  
## Args: see above
## Returns:     list() 
##              X: data.frame with five columns 
##                 ID: subject ID 
##                 x: region-referenced functional predictor
##                 reg: regional argument
##                 func: functional argument
##                 covar: covariate argument
##              y: binary outcomes (vector)
##              lodds: log odds for each subject (vector)
#############################################################################  

# Install missing packages
list.of.packages <- c("MASS", "splines")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) 
    
# Load packages  
library(splines)
library(MASS)

#############################################################################
# 1. Setup evaluation grid, subject-specific score covariance matrix, and basis functions
#############################################################################

g <- expand.grid(w, a) # grid for functional and covariate domain evaluations
S = diag(1 - rho, R) + rho # covariance matrix for subject-specific scores (compound symmetric)
Basisw <- bs(w, df = 5, intercept = TRUE) # functional basis 

#############################################################################
# 2. Simulate subject-specific functional predictors
#############################################################################

X <- NULL # initialize storage for functional predictors

# Simulate functional predictor and non-functional covariate for each subject
for(i in 1:n){ 
  ai <- sample(a, 1) # generate subject-specific covariate 
  temp <- matrix(rep(0, length(w) * R), nrow = length(w), ncol = R) # form temporary storage
  for(k in 1:dim(Basisw)[2]){ # form functional predictor 
    coefs <- mvrnorm(dim(Basisw)[2], rep(0, R), S) # simulate vector of subject-specific scores for kth basis function 
    temp <- temp + Basisw[, k] %*% t(coefs[k, ]) # update functional predictor
  }
  X <- rbind(X, cbind(rep(i, length(temp)), c(temp), 
             rep(1:R, each = length(w)), rep(w, R), rep(ai, length(temp))))
}

X <- data.frame(X)
colnames(X) <- c("ID", "x", "reg", "func", "covar")

#############################################################################
# 3. Define true regression function
#############################################################################

beta <- array(NA, dim = c(R, length(w), length(a)))
  for(r in 1:R){
    if(r %in% c(1:8)){
       beta[r, , ] <- matrix(do.call(mapply, c(function(a, b)  (-1) ^ r * 2 * cos(r * pi /6 * a + pi * b), unname(g))), nrow = length(w))
    } else{
       beta[r, , ] <- matrix(do.call(mapply, c(function(a, b)  (-1) ^ r * 2 * sin((r - R/2) * pi /6 * a + pi * b), unname(g))), nrow = length(w))
    }
}

#############################################################################
# 4. Generate subject-specific probabilities and outcomes
#############################################################################

id <- unique(X$ID) 
by <- w[2] - w[1] # width between functional grid points
y <- rep(NA, n) # storage for outcomes
lodds <- rep(NA, n) # storage for log-odds

# Generate subject-specific probabilities
for(i in 1:n){
  temp <- X[which(X$ID == id[i]), ] 
  betat <- NULL 
  index <- which(a == temp$covar[1]) 
  for(r in 1:R){ # vectorized regression function evaluated at a_i
    betat <- c(betat, beta[r, , index]) 
  }
  lodds[i] <- by * temp$x %*% betat # form logit of subject-specific probabilites 
  y[i] <- rbinom(1, 1, exp(lodds[i])/(1 + exp(lodds[i]))) # generate subject-specific outcomes
}

outcomes <- list(X, y, lodds)
names(outcomes) <- c("X", "y", "lodds")

return(outcomes)

}
