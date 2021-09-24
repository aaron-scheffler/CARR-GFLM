# CARR-GFLM
Code for implementing a covariate-adjusted region-referenced generalized functional linear model.

CONTENTS OF THIS FOLDER
——————————————	

* CARR_GFLM_tutorial.R :           A step-by-step implementation of the covariate-adjusted region-referenced
					generalized functional linear model (CARR-GFLM) and associated procedures described
					in "Covariate-Adjusted Region-Referenced Generalized Functional Linear Model for EEG Data"
					by Scheffler et al. (2019). The procedure assumes that functional predictor is observed on  
					a common grid in the regional, functional, and covariate dimensions.

* CARR_GFLM_data.R : 	      	Function for simulating data for the CARR_GFLM.R function as described
             			           in the simulation section found in Section 4 of "Covariate-Adjusted 
                                                      Region-Referenced Generalized Functional Linear Model for EEG Data" 
                                                      by Scheffler et al. (2019).  

* CARR_GFLM.R : 			Function for fitting CARR-GFLM. Returns the estimated regression function along with Bayesian 
                                                      point-wise confidence intervals and predictions.

* mGFLM.R:                                  Function for fitting a multivariate generalized functional linear model (m-GFLM).

* mGFLMi.R:                                 Function for fitting a m-GFLM with covariate interaction (m-GFLMi).

		 
INTRODUCTION
——————————————	

The contents of this folder allow for implementation of the covariate-adjusted region-referenced generalized functional linear model (CARR-GFLM) 
described in "Covariate-Adjusted Region-Referenced Generalized Functional Linear Model for EEG Data" by Scheffler et al. (2019). Users can simulate 
a sample data frame (CARR_GFLM_data.R) and apply the proposed CARR_GFLM (CARR_GFLM.R). Detailed instructions on how to perform the 
aforementioned procedures and visualize results are included in CARR_GFLM_Tutorial.R. In the tutorial, predictive performance of the proposed model 
is compared to a multivariate generalized functional linear model (m-GFLM) and m-GFLM with covariate interaction (m-GFLMi).

REQUIREMENTS
——————————————	

The included R programs require R 3.5.1 (R Core Team, 2018) and the packages listed in CARR_GFLM_Tutorial.R.


INSTALLATION
——————————————	

Load the R program files into the global environment and install required packages using commands in CARR_GFLM_Tutorial.R.
