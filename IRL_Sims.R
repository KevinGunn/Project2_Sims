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
# $ git pull origin master
# Pushes the changes in your local repository up to the remote repository 
# you specified as the origin

# Libraries
library(MASS)
library(plotly)
library(dfoptim)
library(quadprog)

#set seed
set.seed(12345)


# Penalty Function and linear decision rule functions:

norm <- function(x) sqrt(sum(x^2))

hinge <- function(s){ ifelse(s>0,s,0) } 

trt_rule <- function(z){ ifelse(z>0,1,0) }

empirical_value <- function(rule,trt,prop,outcome){
  
  mean(ifelse(rule==trt,1,0)*prop^(-1)*outcome)
  
}

# Obtain eta^{opt} and assume it is the clinicians decision rule.
clin_eta_opt <- function( eta, X, params_list, M, lambda){
  
  etax <- as.vector(eta %*% t(X))
  
  EY <- mean( params_list[[5]] + as.vector( params_list[[1]] %*% t(X) ) + 
                trt_rule(etax) * as.vector( params_list[[2]] %*% t(X) ) ) 
  
  EZ_lambda <-  mean( params_list[[6]] + as.vector( params_list[[3]] %*% t(X) ) + 
                        trt_rule(etax) * as.vector( params_list[[4]] %*% t(X) ) ) - lambda
  
  # Penalized Value function.
  V_pen <- EY + M * hinge(EZ_lambda)
  
  #print(EY);print(hinge(EZ_lambda));print(V_pen)
  
  return(V_pen)
}

# Value function estimation
VY_est <- function(rule){ mean( QY$coefficients[1] + 
                                  as.vector(QY$coefficients[2:(1+p)] %*% t(X)) +
                                  rule*as.vector(QY$coefficients[-(1:(1+p))] %*% t(X)) ) 
}

VZ_tol_est <- function(rule, tol){
  
  VZ = mean( QZ$coefficients[1] + as.vector(QZ$coefficients[2:(1+p)] %*% t(X)) + rule*as.vector(QZ$coefficients[-(1:(1+p))] %*% t(X)) )
  return(VZ - tol)                            
}

# Obtain value function for eta^{opt} estimate.
value_eta <- function( eta, X, params_list, M, lambda){
  
  etax <- as.vector(eta %*% t(X))
  
  EY <- mean( params_list[[5]] + as.vector( params_list[[1]] %*% t(X) ) + 
                trt_rule(etax) * as.vector( params_list[[2]] %*% t(X) ) ) 
  
  EZ_lambda <-  mean( params_list[[6]] + as.vector( params_list[[3]] %*% t(X) ) + 
                        trt_rule(etax) * as.vector( params_list[[4]] %*% t(X) ) ) - lambda
  
  # Penalized Value function.
  V_pen <- EY + M * hinge(EZ_lambda)
  
  return(V_pen)
}


# Provide Settings:

p <- 5
n <- 1000
N <- 10000
sims <- 500

#params

theta_0Y <- c(1,0,0,0,1)
theta_1Y <- c(0.25, 0.25, 0.25, 0, 0)
theta_0Z <- c(1,0,1,0,0)
theta_1Z <- c(1,0,0,2,1)

params <- list(theta_0Y, theta_1Y, theta_0Z, theta_1Z,0,0)

# Generate Covariates Data 
X <- mvrnorm(n=N, mu=rep(0,p), diag(p) )

#######################
# Generate a reasonable value of lambda
# Get Contrast function
# C_X = as.vector( params_list[[4]] %*% t(X) )  - as.vector( params_list[[3]] %*% t(X) )
# VZ =  mean( as.vector( params_list[[3]] %*% t(X) ) + 
#             trt_rule(C_X) * as.vector( params_list[[4]] %*% t(X) ) )
#######################

########################################################################

## Plot all combinations of M and lambda to get a reasonable value for their optimal values.

M <- seq(0.5,5,0.5)
lambda <- seq(0.1,3,0.1)

grid <- expand.grid(M,lambda)
V <- rep(0,nrow(grid))
eta_opt_mat <- matrix(0,ncol=5,nrow=nrow(grid))
EY_EZ <- matrix(0,ncol=2,nrow=nrow(grid))

for(i in 1:nrow(grid)){
  
  eta_func <- function(eta){ clin_eta_opt(eta=eta, X=X, params_list = params,
                                          M=grid[i,1], lambda = grid[i,2])}
  
  eta_initial <- rep(0,5)
  

  eta_opt <- nmk(par = eta_initial, 
                 eta_func,
                 control=list(maximize=TRUE)
  )
  
  print(eta_opt$message);print(eta_func(eta_opt$par))
  eta_opt_mat[i,] <- eta_opt$par
  V[i] <- eta_func(eta_opt$par)
  EY_EZ[i,1] <- mean( params[[5]] + as.vector( params[[1]] %*% t(X) ) + 
                        trt_rule(as.vector(eta_opt$par %*% t(X))) * 
                        as.vector( params[[2]] %*% t(X) ) ) 
  
  EY_EZ[i,2] <-  hinge(mean( params[[6]] + as.vector( params[[3]] %*% t(X) ) + 
                               trt_rule(as.vector(eta_opt$par %*% t(X))) * 
                               as.vector( params[[4]] %*% t(X) ) ) - grid[i,2])
  
  
}

grid_vals <- cbind(grid,V); colnames(grid_vals) <- c("M","lambda","V")
gv_df <- as.data.frame(grid_vals)

# Scatterplot
plot_3d <- plot_ly(gv_df, x = ~M, y = ~lambda, z = ~ V, 
                   marker = list(color = ~ V, colorscale = c('#FFE1A1', '#683531'), 
                                 showscale = TRUE)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'M'),
                      yaxis = list(title = 'Lambda'),
                      zaxis = list(title = 'V')),
         
         annotations = list(
           x = 1,
           y = 0.5,
           text = 'Value Function',
           xref = 'paper',
           yref = 'paper',
           showarrow = FALSE
         ))

plot_3d


# Based on scatterplot let M=3.7, lambda=0.2, and V = 0.1164.
EY_EZ
M_opt = gv_df[167,1]; lambda_opt = gv_df[167,2]; eta_opt <- eta_opt_mat[167,]
#M_opt = gv_df[3,1]; lambda_opt = gv_df[3,2]; eta_opt <- eta_opt_mat[3,]

#################################################################################

#Generation of Data assuming clinician makes optimal decision.

# Generate Covariates Data 
X <- mvrnorm(n=n, mu=rep(0,p), diag(p) )

# 100% Correct decisions
opt_model <- as.vector(eta_opt %*% t(X))

clin_model <- opt_model

################
# Add error to optimal model so clinician is slightly imperfect.
#clin_model <- opt_model + rnorm(n,0,1.5)

# Clinician is approximately 91% correct.
#sum(trt_rule(opt_model) == trt_rule(clin_model))  
###############

# Treatment
A <- trt_rule(clin_model)

# Main Outcome of Interest
Y <- as.vector( theta_0Y %*% t(X) ) + A * as.vector( theta_1Y %*% t(X) ) + rnorm(n,sd=0.5)

# Risk Outcome
Z <- as.vector( theta_0Z %*% t(X) ) + A * as.vector( theta_1Z %*% t(X) ) + rnorm(n,sd=0.5)


###################################################################################



# Q-model to find theta values

QY <- lm(Y ~ X + A:X)
QZ <- lm(Z ~ X + A:X)

params_est <- list(QY$coefficients[2:(p+1)], QY$coefficients[-(1:(p+1))], 
                   QZ$coefficients[2:(p+1)], QZ$coefficients[-(1:(p+1))],
                   QY$coefficients[1],QZ$coefficients[1])


# Start algorithm.
# Create empty matrices to store values in.
lin_rules_mat <- matrix(0,ncol=5,nrow=sims)


#Step 1 - Grid Search
#lambda_est <- c(1,1.7,2.0)
lambda_est <- c(0.15,0.25,0.05)


eta_mat <- matrix(0,ncol=p,nrow=50)
k_eta_mat_lambda <- matrix(0, ncol=p, nrow = length(lambda_est))

M_est <- rep(0, length(lambda_est) )
V_est <- rep(0, length(lambda_est) )


for (i in 1:length(lambda_est)) {
  
  # Clinicians Value function
  V_Y_clin <-  mean(Y)
  V_Z_clin_tol <- mean(Z) - lambda_est[i]  #EY_EZ[167,2] - lambda_est[i]
  
  mu_clin <- c( V_Y_clin, hinge(V_Z_clin_tol) )
  
  # Initialize decision rule
  lin_rule_coeff <- c(0,0,0,0,0)
  lin_mod1 <- as.vector(lin_rule_coeff %*% t(X))
  rule1 <- trt_rule(lin_mod1)
  
  VY_learner <- VY_est(rule1)
  VZ_learner <- VZ_tol_est(rule1, lambda_est[i])
  
  mu_learner <- c(VY_learner, hinge(VZ_learner) )
  
  # Begin AL-IRL
  k = 1;q=1;eps = 0.0001
  
  eta_opt_k <- lin_rule_coeff
  
  # Create data matrix for svm.
  mu_mat <- rbind(mu_clin, mu_learner)
  colnames(mu_mat) <- c("VY","VZ_tol")
  rownames(mu_mat) <- NULL
  
  Labels <- c(1,-1)
  L_mu_mat <- cbind(Labels,mu_mat)
  
  while(k <= 50){
    
    # QP set up.
    m = length(mu_clin)
    Dm = diag(m+1) # min ||w||^2
    Dm[nrow(Dm)-2, ncol(Dm)-2] <- Dm[nrow(Dm)-1, ncol(Dm)-1] <- 0.0000001
    dv = rep(0,m+1)
    
    c_m1 <- c(1,1,0)
    #c_m2 <- rbind(c(1,0,0),c(0,1,0))
    
    VY_diff <- L_mu_mat[1,2] - L_mu_mat[-1,2]
    VZ_diff <- L_mu_mat[1,3] - L_mu_mat[-1,3]
    V_diff <- cbind(VY_diff,VZ_diff,-q)
    
    
    Am <- rbind(c_m1, V_diff) #,c_m2 )
    
    bv <- c(1,rep(0,k))
    sol <- solve.QP(Dm,dv,t(Am),bv,meq=1)

    Mk = sol$solution[2]/sol$solution[1]
    qk = sol$solution[3]
     
    
    # Plot of Value functions.
    plot(L_mu_mat[,-1],col="blue", pch=19, main="Value Function Estimation",cex=1.2)
    points(x= L_mu_mat[1,2],y= L_mu_mat[1,3],bg="red",pch=22,cex=1.5) 
    # Hyperplane
    #abline( b=1/Mk , col="green", lty=4)
    
  
    # Obtain estimate of eta^opt
    eta_func <- function(eta){ value_eta(eta, X=X, params_list = params_est,
                                         M=Mk, lambda = lambda_est[i]) }
    
    eta_opt_k <- nmk(par = rep(0,5) , eta_func, control=list(maximize=TRUE))$par
    
    # Get decision rule.
    lin_modk <- as.vector(eta_opt_k %*% t(X))
    rulek <- trt_rule(lin_modk)
    
    VY_learner <- VY_est(rulek)
    VZ_learner <- VZ_tol_est(rulek, lambda_est[i])
    
    # Update data
    mu_learner <- c(VY_learner, hinge(VZ_learner) )
    L_mu_mat <- rbind(L_mu_mat, c(-1,mu_learner) )
    
    
    #sum(rulek == A)
    # update iteration step.
    if(k > 5 & abs(qk)<eps){break}
    
    k = k+1
    
  }
  
  M_est[i] <- Mk
  k_eta_mat_lambda[i,] <-  eta_opt_k
  V_est[i] <-  mu_learner[1] + Mk*mu_learner[2]
  
}




