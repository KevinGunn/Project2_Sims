###################################################
## Kevin Gunn                                    ##
## 4/17/2019                                     ##
## Project 2 - Constrained IRL                   ##
##                                               ##
###################################################

# Git: https://github.com/KevinGunn/Project2_Sims.git
# Git Instructions:

#Step 1 - Update local repository

#$ cd /path/to/my/codebase
#$ git init      (1)
#$ git add .     (2)
#$ git commit -m "comments"    (3)

# Step 2 - Push to online repository.

# git remote add origin https://github.com/KevinGunn/Project2_Sims.git (1)
# $ git push origin master (2)
# Pushes the changes in your local repository up to the remote repository 
# you specified as the origin


# Libraries

library(MASS)

#set seed
set.seed(12345)


# Penalty Function and linear decision rule functions:

hinge <- function(s){ ifelse(s>0,s,0) } 

trt_rule <- function(z){ ifelse(z>0,1,0) }

# Provide Settings:

p <- 5
n <- 1000

theta_0Y <- c(1,0,0,0,1)
theta_1Y <- c(1, 1, 1, 0, 0)
theta_0Z <- c(0.25,0.25,0,0,0.25)
theta_1Z <- rep(0.5,5)

#params
params <- list(theta_0Y, theta_1Y, theta_0Z, theta_1Z)

# Generate Covariates Data 
X <- mvrnorm(n=n, mu=rep(0,p), diag(p) )

#######################
# Generate a reasonable value of lambda
# Get Contrast function
# C_X = as.vector( params_list[[4]] %*% t(X) )  - as.vector( params_list[[3]] %*% t(X) )
# VZ =  mean( as.vector( params_list[[3]] %*% t(X) ) + 
#             trt_rule(C_X) * as.vector( params_list[[4]] %*% t(X) ) )
#######################

M <- 1.5
lambda <- 0.3

# Obtain eta^{opt} and assume it is the clinicians decision rule.
clin_eta_opt <- function( eta, X, params_list, M, lambda){
  
  etax <- as.vector(eta %*% t(X))
  
  EY <- mean( as.vector( params_list[[1]] %*% t(X) ) + 
              trt_rule(etax) * as.vector( params_list[[2]] %*% t(X) ) ) 
  
  EZ_lambda <-  mean( as.vector( params_list[[3]] %*% t(X) ) + 
                trt_rule(etax) * as.vector( params_list[[4]] %*% t(X) ) ) - lambda

  # Penalized Value function.
  V_pen <- EY - M * hinge(EZ_lambda)
  
  #print(EY);print(hinge(EZ_lambda))
  
  return(V_pen)
}


eta_func <- function(eta){clin_eta_opt(eta=eta, X=X, params_list = params,
                                      M=M, lambda = lambda)}

eta_opt <- optim(par = c(1,1,1,0,0), 
                  eta_func,
                  method = "BFGS"
                )$par

patient_rule <- as.vector(eta_opt %*% t(X))


# Treatment
A <- trt_rule(patient_rule)

# Main Outcome of Interest
Y <- as.vector( theta_0Y %*% t(X) ) + A * as.vector( theta_1Y %*% t(X) ) + rnorm(n)

# Risk Outcome
Z <- as.vector( theta_0Z %*% t(X) ) + A * as.vector( theta_1Z %*% t(X) ) + rnorm(n)





