###############################################################
### Sample code for generating and estimating               ###
### simulated datasets.                                     ###
###############################################################
### The list of variables used in this analysis.
### n: Number of subjects in each cluster.
### K: The number of variables.
### K1: The number of signal variables.
### prior_mu: The prior information of variables.
### prior_sigma: The prior information of interactions.


###############################################################
### Call packages and functions:
library(MASS)
library(dplyr)
library(sparcl)
library(cluster)
library(factoextra)
library(CEoptim)
library(mclust)
library(mvtnorm)
library(JGL)
library(glasso)
library(Matrix)

source("Functions.R")


###############################################################
### Generate simulated datasets:    

set.seed(11)
simul_data <- PECM_simu_gene_2(K = 1000,K1 = 20)

# Set values for input variables.
X <- simul_data$X
prior_mu <- c(1:10,11:20)
prior_sigma <- matrix(c(61:79,62:80,62:80,61:79),ncol = 2)

# First fully trust the prior information
aexci <- APECM(X,2,prior_mu,prior_sigma,3,0.5,5)
prior_information <- aexci$label
# Calculate the CER
cal_CER(prior_information,simul_data$label)

# Then use the proposed method
excit <- PECM(X,2,prior_information = prior_information,lambda1 = 5,lambda2 = 0.4,eta = 0.1,B = 10)
# Calculate the CER
cal_CER(excit$label,simul_data$label)

