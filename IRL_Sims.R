###################################################
## Kevin Gunn                                    ##
## 4/17/2019                                     ##
## Project 2 - Constrained IRL                   ##
##                                               ##
###################################################

# Git: https://github.com/KevinGunn/Project2_Sims.git
# Git Instructions:

#Step 1 - Make/Update local repository

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
library(plotly)
library(e1071)

#set seed
set.seed(12345)


# Penalty Function and linear decision rule functions:

hinge <- function(s){ ifelse(s>0,s,0) } 

trt_rule <- function(z){ ifelse(z>0,1,0) }

empirical_value <- function(rule,trt,prop,outcome){
  
  mean(ifelse(rule==trt,1,0)*prop^(-1)*outcome)
  
}

# Provide Settings:

p <- 5
n <- 1000
sims <- 500

#params

theta_0Y <- c(1,0,0,0,1)
theta_1Y <- c(0.2, 0.25, 0.25, 0, 0)
theta_0Z <- c(1,0.25,1,0,0.25)
theta_1Z <- c(0.5,0.5,0.5,-0.5,1)

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

########################################################################

## Plot all combinations of M and lambda to get a reasonable value for their optimal values.

M <- seq(0.1,5,0.1)
lambda <- seq(0.1,3,0.1)

grid <- expand.grid(M,lambda)
V <- rep(0,nrow(grid))
eta_opt_mat <- matrix(0,ncol=5,nrow=nrow(grid))

for(i in 1:nrow(grid)){
    
  # Obtain eta^{opt} and assume it is the clinicians decision rule.
  clin_eta_opt <- function( eta, X, params_list, M, lambda){
    
    etax <- as.vector(eta %*% t(X))
    
    EY <- mean( as.vector( params_list[[1]] %*% t(X) ) + 
                trt_rule(etax) * as.vector( params_list[[2]] %*% t(X) ) ) 
    
    EZ_lambda <-  mean( as.vector( params_list[[3]] %*% t(X) ) + 
                  trt_rule(etax) * as.vector( params_list[[4]] %*% t(X) ) ) - lambda
  
    # Penalized Value function.
    V_pen <- EY - M * hinge(EZ_lambda)
    
    #print(EY);print(hinge(EZ_lambda));print(V_pen)
    
    return(V_pen)
  }
  
  
  eta_func <- function(eta){ clin_eta_opt(eta=eta, X=X, params_list = params,
                                        M=grid[i,1], lambda = grid[i,2])}
  
  eta_initial <- rep(0,5)
  
  eta_opt <- optim(par = eta_initial, 
                    eta_func,
                    control=list(fnscale=-1),
                    method = "BFGS"
                  )
  
  #print(eta_opt);print(eta_func(eta_opt$par))
  eta_opt_mat[i,] <- eta_opt$par
  V[i] <- eta_func(eta_opt$par)
  
}

grid_vals <- cbind(grid,V); colnames(grid_vals) <- c("M","lambda","V")
gv_df <- as.data.frame(grid_vals)

# Scatterplot
plot_3d <- plot_ly(gv_df, x = ~M, y = ~lambda, z = ~-V, 
             marker = list(color = ~ V, colorscale = c('#FFE1A1', '#683531'), 
                           showscale = TRUE)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'M'),
                      yaxis = list(title = 'Lambda'),
                      zaxis = list(title = '-V')),
  
    annotations = list(
          x = 1,
          y = 0.5,
          text = 'Negative of Value Function',
          xref = 'paper',
          yref = 'paper',
          showarrow = FALSE
    ))

plot_3d


# Based on scatterplot let M=3.7, lambda=0.2, and V = 0.1164.
M_opt = 3.7; lambda_opt = 0.2; eta_opt <- eta_opt_mat[87,]

#################################################################################

#Generation of Data assuming clinician makes optimal decision.

# 100% Correct decisions
opt_model <- as.vector(eta_opt %*% t(X))

# Add error to optimal model so clinician is slightly imperfect.
clin_model <- opt_model + rnorm(n,0,1.5)

# Clinician is approximately 91% correct.
sum(trt_rule(opt_model) == trt_rule(clin_model))  

# Treatment
A <- trt_rule(clin_model)

# Main Outcome of Interest
Y <- as.vector( theta_0Y %*% t(X) ) + A * as.vector( theta_1Y %*% t(X) ) + rnorm(n)

# Risk Outcome
Z <- as.vector( theta_0Z %*% t(X) ) + A * as.vector( theta_1Z %*% t(X) ) + rnorm(n)


###################################################################################

# Fit propensity score
prop_score_model <- glm(A ~ X, family="binomial")

prop_score <- prop_score_model$fitted.values

# Q-model to find theta values

QY <- lm(Y ~ X + A:X)
QZ <- lm(Z ~ X + A:X)

VY_est <- function(rule){ mean( QY$coefficients[1] + 
                      as.vector(QY$coefficients[2:(1+p)] %*% t(X)) +
                      rule*as.vector(QY$coefficients[-(1:(1+p))] %*% t(X)) ) 
}

VZ_tol_est <- function(rule, tol){
  
  VZ = mean( QZ$coefficients[1] + as.vector(QZ$coefficients[2:(1+p)] %*% t(X)) + rule*as.vector(QZ$coefficients[-(1:(1+p))] %*% t(X)) )
  return(VZ - tol)                            
}

# Start algorithm.
# Create empty matrices to store values in.

lin_rules_mat <- matrix(0,ncol=5,nrow=sims)

#Step 1 - Grid Search
lambda_est <- seq(0.1,0.3,0.05)


for (i in 1:length(lambda_est)) {
  
  # Clinicians Value function
  V_Y_clin <- mean(Y)
  V_Z_clin_tol <- mean(Z) - lambda_est[i]

  mu_clin <- c( V_Y_clin, - hinge(V_Z_clin_tol) )
  
  # Initialize decision rule
  lin_rule_coeff <- rep(0.01,5)
  lin_mod1 <- as.vector(lin_rule_coeff %*% t(X))
  rule1 <- trt_rule(lin_mod1)
  
  VY_learner <- VY_est(rule1)
  VZ_learner <- VZ_tol_est(rule1, lambda_est[i])
  
  mu_learner <- c(VY_learner, - hinge(VZ_learner) )
  
  # Begin AL-IRL
  k = 1;q=1;eps = 0.01
  
  # Create data matrix for svm.
  mu_mat <- rbind(mu_clin, mu_learner)
  Labels <- c(1,-1)
  L_mu_mat <- cbind(Labels,mu_mat)
  
  while(q > eps){
                
    classifier = svm(formula = Purchased ~ ., 
                 data = training_set, 
                 type = 'C-classification', 
                 kernel = 'linear') 

  }
  
}
