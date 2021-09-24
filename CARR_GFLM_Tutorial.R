## CARR_GFLM_Tutorial.R
################################################################################################
## Description: A step-by-step implementation of the covariate-adjusted region-referenced
## generalized functional linear model (CARR-GFLM) and associated procedures described
## in "Covariate-Adjusted Region-Referenced Generalized Functional Linear Model for EEG Data"
## by Scheffler et al. (2019). The procedure assumes that functional predictor is observed on  
## a common grid in the regional, functional, and covariate dimensions. Predictive performance of
## the proposed model is compared to a multivariate generalized functional linear model (m-GFLM) 
## and m-GFLM with covariate interaction (m-GFLMi).
################################################################################################
## Functions implemented: 
## CARR_GFLM.R, CARR_GFLM_data.R, mGFLM.R, mGFLMi.R
################################################################################################
## Tutorial Outline:
##    1. Form region-referenced functional predictors and simulate binary outcomes (CARR_GFLM_data.R)
##    2. Fit CARR-GFLM, m-GFLM, and m-GFLMi to simulated data (CARR_GFLM.R, mGFLM.R, mGFLMi.R)
##    3. Visualize the true and estimated regression function
##    4. Compare predictive performance between CARR-GFLM and m-GFLM, m-GFLMi
################################################################################################

# Install missing packages
list.of.packages <- c("mgcv", "splines", "data.table", "ppls", "pracma", "ggplot2", "latex2exp", "plotROC")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) 
    
# Load packages  
library(data.table)
library(splines)
library(mgcv)
library(ppls)
library(pracma)
library(ggplot2)
library(latex2exp)
library(plotROC)


################################################################################################
# 1. Form region-referenced functional predictors and simulate binary outcomes (CARR_GFLM_data.R)
################################################################################################

# Domains
w <- 1:50/50 # functional grid points
a <- 1:30/30 # covariate grid points
R <- 15 # number of regions

# Simulate data from simulation design at N = 500 and \rho = .1
simDat <- CARR_GFLM_data(n = 500,                 # sample size
                         rho = .1,                # level of regional dependency
                         R = R,                   # number of regions
                         w = w,                   # functional domain
                         a = a)                   # covariate domain

################################################################################################
# 2. Fit CARR-GFLM to simulated data (CARR_GFLM.R)
################################################################################################

# Marginal basis functions
Basisw <- bs(w, df = 5, intercept = TRUE) # functional basis 
Basisw <- suppressWarnings(apply(as.matrix(Basisw), 2, function(x) x/sqrt(x %*% x))) # normalize

Basisa <- bs(a, df = 5, intercept = TRUE) # covariate basis 
Basisa <- suppressWarnings(apply(as.matrix(Basisa), 2, function(x) x/sqrt(t(x) %*% x))) # normalize

Basisr <- cbind(rep(1, R), eye(R)) # regional basis 
Basisr <- suppressWarnings(matrix(apply(as.matrix(Basisr), 2, function(x) x/sqrt(t(x) %*% x)), ncol = dim(Basisr)[2])) # normalize

# Marginal penalty matrices
Pw <- Penalty.matrix(5, order = 1) # functional penalty matrix
Pa <- Penalty.matrix(5, order = 1) # covariate penalty matrix
Pr <- diag(c(0, rep(1, R))) # regional penalty matrix

# Fit CARR-GFLM 
carr_gflm <- CARR_GFLM(X = simDat$X, # functional predictors
                  y = simDat$y,      # outcome 
                  w = w,             # functional domain 
                  a = a,             # covariate domain 
                  Basisw = Basisw,   # functional basis 
                  Pw = Pw,           # functional penalty
                  Basisa = Basisa,   # covariate basis 
                  Pa = Pa,           # covariate penalty 
                  Basisr = Basisr,   # regional basis
                  Pr = Pr,           # regional penalty 
                  dist = "binomial", # distribution of outcome
                  pvar = .95,        # proportion of variation 
                  ci = 1,            # return confidence intervals (1 = yes, 0 = no)
                  kbasis = c(10, 5), # degrees of freedom for region-referenced mean surface smoothing, c(functional, covariate)
                  weights = 0        # implement weighted smoothing for region-referenced mean surface (1 = yes, 0 = no), for binomial outcomes only
                 )

# Fit m-GFLM
m_gflm <- mGFLM  (X = simDat$X,      # functional predictors
                  y = simDat$y,      # outcome 
                  w = w,             # functional domain 
                  Basisw = Basisw,   # functional basis 
                  Pw = Pw,           # functional penalty
                  dist = "binomial", # distribution of outcome
                  pvar = .95,        # proportion of variation 
                  kbasis = c(10),    # degrees of freedom for region-referenced mean surface smoothing, c(functional, covariate)
                  weights = 0        # implement weighted smoothing for region-referenced mean surface (1 = yes, 0 = no), for binomial outcomes only
                 )

# Fit m-GFLMi
m_gflmi <- mGFLMi(X = simDat$X,      # functional predictors
                  y = simDat$y,      # outcome 
                  w = w,             # functional domain 
                  Basisw = Basisw,   # functional basis 
                  Pw = Pw,           # functional penalty
                  dist = "binomial", # distribution of outcome
                  pvar = .95,        # proportion of variation 
                  kbasis = c(10),    # degrees of freedom for region-referenced mean surface smoothing, c(functional, covariate)
                  weights = 0        # implement weighted smoothing for region-referenced mean surface (1 = yes, 0 = no), for binomial outcomes only
                 )

################################################################################################
# 3. Visualize the true and estimated regression function
################################################################################################

# True regression function
g <- expand.grid(w, a) 
beta <- array(NA, dim = c(R, length(w), length(a)))
  for(r in 1:R){
    if(r %in% c(1:8)){
       beta[r, , ] <- matrix(do.call(mapply, c(function(a, b)  (-1) ^ r * 2 * cos(r * pi /6 * a + pi * b), unname(g))), nrow = 50)
    } else{
       beta[r, , ] <- matrix(do.call(mapply, c(function(a, b)  (-1) ^ r * 2 * sin((r - R/2) * pi /6 * a + pi * b), unname(g))), nrow = 50)
    }
}

# Plot true (left column) and estimated (right column) regression functions for three regions
par(mfrow=c(3, 2), oma = c(2, 1, 1, 1), mar = c(2, 1, 1, 1))
for(r in c(1, 2, 3)){
  persp(beta[r, , seq(1, 30, by = 3)], theta = 300, main = paste("True: R=", r, sep = ""), 
        zlim = c(- 3, 3), ticktype = "detailed", xlab = "omega", ylab = "a",  zlab = "beta")
  persp(carr_gflm$betaHat[[r]][, seq(1, 30, by = 2)], theta = 300, main = paste("Est: R=", r, sep = ""), 
        zlim = c(-3, 3), ticktype = "detailed", xlab = "omega", ylab = "a",  zlab = "beta")
}

################################################################################################
# 4. Compare predictive performance between CARR-GFLM and m-GFLM, m-GFLMi
################################################################################################

# Form matrix with true outcomes and predicted probabilities 
AUC <- data.frame(cbind(rep(simDat$y, 3), 
                  c(carr_gflm$fit, m_gflm$fit, m_gflmi$fit),
                  rep(c("CARR-GFLM", "m-GFLM", "m-GFLMi"), each = length(simDat$y))))
names(AUC) <- c("true", "fit", "model")   
AUC$fit <- as.numeric(as.character(AUC$fit))
AUC$true <- as.numeric(AUC$true) - 1

# Plot ROC curve
g <- ggplot(data  = AUC, aes(d = true, m = fit, color = model)) + geom_roc(labels = FALSE, n.cuts = 0) +
     style_roc(xlab = "1 - Specificity", ylab = "Sensitivity") +
     ggtitle("ROC Curve by Model")
print(g)


