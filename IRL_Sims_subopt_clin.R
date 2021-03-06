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

# Condition handling code to avoid bad starting values in nmk.
show_condition <- function(code) {
  tryCatch(code,
           error = function(c) "error",
           warning = function(c) "warning",
           message = function(c) "message"
  )
}


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
  
  EZ_lambda <-  mean( params_list[[6]] + as.vector( params_list[[3]] %*% t(X) ) - 
                        trt_rule(etax) * as.vector( params_list[[4]] %*% t(X) ) ) - lambda
  
  # Penalized Value function.
  V_pen <- EY - M * hinge(EZ_lambda)
  
  #print(EY);print(hinge(EZ_lambda));print(V_pen)
  
  return(V_pen)
}

# Value function estimation
VY_est <- function(rule){ mean( QY$coefficients[1] + 
                                  as.vector(QY$coefficients[2:(1+p)] %*% t(X)) +
                                  rule*as.vector(QY$coefficients[-(1:(1+p))] %*% t(X)) ) 
}

VZ_tol_est <- function(rule, tol){
  
  VZ = mean( QZ$coefficients[1] + as.vector(QZ$coefficients[2:(1+p)] %*% t(X)) + 
               rule*as.vector(QZ$coefficients[-(1:(1+p))] %*% t(X)) )
  VZ_tol <- VZ - tol
  
  return(VZ_tol)
}

# Obtain value function for eta^{opt} estimate.
value_eta <- function( eta, X, params_list, M, lambda){
  
  etax <- as.vector(eta %*% t(X))
  
  EY <- mean( params_list[[5]] + as.vector( params_list[[1]] %*% t(X) ) + 
                trt_rule(etax) * as.vector( params_list[[2]] %*% t(X) ) ) 
  
  EZ_lambda <-  mean( params_list[[6]] + as.vector( params_list[[3]] %*% t(X) ) + 
                        trt_rule(etax) * as.vector( params_list[[4]] %*% t(X) ) ) - lambda
  
  # Penalized Value function.
  V_pen <- EY - M * hinge(EZ_lambda)
  
  return(V_pen)
}


# Provide Settings:

p <- 5
n <- 2000
N <- 100000

#params

theta_0Y <- c(1,0,1,0,1)
theta_1Y <- c(1, 0.5, 0.5, 1, 1)
theta_0Z <- c(0,0,1,1,0)
theta_1Z <- c(0.25,1,1,0,1)

params <- list(theta_0Y, theta_1Y, theta_0Z, theta_1Z,0,0)

# Generate Covariates Data 
X <- mvrnorm(n=N, mu=rep(0,p), diag(p) )

####################################################################################
# Generate a reasonable value of lambda

opt_rule1 <- as.vector(params[[2]] %*% t(X))

opt_rule2 <- as.vector(params[[4]] %*% t(X))

sum(trt_rule(opt_rule1) == trt_rule(opt_rule2)) / length(opt_rule1)


EY <- mean( params[[5]] + as.vector( params[[1]] %*% t(X) ) + 
              trt_rule(opt_rule1) * as.vector( params[[2]] %*% t(X) ) ) 

EZ_optr1 <- mean( params[[6]] + as.vector( params[[3]] %*% t(X) ) - 
                    trt_rule(opt_rule1) * as.vector( params[[4]] %*% t(X) ) )

EZ_optr2 <-  mean( params[[6]] + as.vector( params[[3]] %*% t(X) ) - 
                     trt_rule(opt_rule2) * as.vector( params[[4]] %*% t(X) ) )

# Expected Value Estimates.
EY
EZ_optr1
EZ_optr2

lambda_opt <- -0.75

EZ_optr2 - lambda_opt
EZ_optr1 - lambda_opt
#####################################################################################

#####################################################################################

## Plot all combinations of M and lambda to get a reasonable value for their optimal values.

M <- seq(0.1,1,0.1)
lambda <- lambda_opt #seq(0.1,3,0.1)

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
                        trt_rule(as.vector(eta_opt$par %*% t(X))) * as.vector( params[[2]] %*% t(X) ) ) 
  
  EY_EZ[i,2] <-  hinge(mean( params[[6]] + as.vector( params[[3]] %*% t(X) ) - 
                               trt_rule(as.vector(eta_opt$par %*% t(X))) * as.vector( params[[4]] %*% t(X) ) ) - grid[i,2])
  
  
}

grid_vals <- cbind(grid,V); colnames(grid_vals) <- c("M","lambda","V")
gv_df <- as.data.frame(grid_vals)


# Based on scatterplot let M=0.8
EY_EZ
eta_opt_mat
M_opt = gv_df[8,1]; lambda_opt = gv_df[10,2]; eta_opt <- eta_opt_mat[8,]

#################################################################################

opt_rule_eta <- as.vector(eta_opt %*% t(X))

EZ_opt <-  mean( params[[6]] + as.vector( params[[3]] %*% t(X) ) - 
                   trt_rule(opt_rule_eta) * as.vector( params[[4]] %*% t(X) ) )

# Expected Value Estimates.
EY
EZ_opt


#################################################################################

# number of simulations is set for 500.
num_sims <- 500
lambda_est <- lambda_opt

#Store output in following matrices/vectors.
eta_mat_lambda <- matrix(0, ncol=p, nrow = num_sims)

M_est_vec <- rep(0, num_sims )
V_est_vec <- rep(0, num_sims )
V_clin_vec <- rep(0, num_sims )

VY_est_vec <- rep(0, num_sims )
VZ_est_vec <- rep(0, num_sims )


VY_clin_vec <- rep(0, num_sims )
VZ_clin_vec <- rep(0, num_sims )

clin_v_all <- matrix(0,ncol=2,nrow=num_sims)

PCD_vec <- rep(0, num_sims )
A_vec <- rep(0, num_sims )

for(sim in 1:num_sims){
  
  # Generate Covariates Data 
  X <- mvrnorm(n=n, mu=rep(0,p), diag(p) )
  
  # 100% Correct decisions
  opt_model <- as.vector(eta_opt %*% t(X))
  
  #suboptimal clinical decision
  eta_sub_opt <- eta_opt + rnorm(5 , 0 , 0.25 )
  clin_model <- as.vector(eta_sub_opt %*% t(X))
  
  # Treatment
  A <- trt_rule(clin_model)
  A_true <- trt_rule(opt_model)
  
  # Main Outcome of Interest
  Y <- as.vector( theta_0Y %*% t(X) ) + A * as.vector( theta_1Y %*% t(X) ) + rnorm(n,sd=0.5)
  
  # Risk Outcome
  Z <- as.vector( theta_0Z %*% t(X) ) - A * as.vector( theta_1Z %*% t(X) ) + rnorm(n,sd=0.5)
  
  
  ###################################################################################
  # Q-model to find theta values
  
  QY <- lm(Y ~ X + A:X)
  QZ <- lm(Z ~ X + A:X)
  
  params_est <- list(QY$coefficients[2:(p+1)], QY$coefficients[-(1:(p+1))], 
                     QZ$coefficients[2:(p+1)], QZ$coefficients[-(1:(p+1))],
                     QY$coefficients[1],QZ$coefficients[1])
  
  
  # Clinicians Value function
  V_Y_clin <-  mean(Y)
  V_Z_clin_tol <- mean(Z) - lambda_est
  
  mu_clin <- c( V_Y_clin, -hinge(V_Z_clin_tol) )
  
  clin_v_all[sim,] <- mu_clin
  
  # Initialize decision rule
  lin_rule_coeff <- c(1,0,0,1,1)
  lin_mod1 <- as.vector(lin_rule_coeff %*% t(X))
  rule1 <- trt_rule(lin_mod1)
  
  VY_learner <- VY_est(rule1)
  VZ_learner <- VZ_tol_est(rule1, lambda_est)
  
  mu_learner <- c(VY_learner, -hinge(VZ_learner) )
  
  # Begin AL-IRL Algorithm.
  k = 1;eps = 0.00001
  
  eta_opt_k <- lin_rule_coeff
  
  # Create data matrix for svm.
  mu_mat <- rbind(mu_clin, mu_learner)
  colnames(mu_mat) <- c("VY","VZ_tol")
  rownames(mu_mat) <- NULL
  
  Labels <- c(1,-1)
  L_mu_mat <- cbind(Labels,mu_mat)
  
  while(k <= 20){
    
    q = rep(1,k)
    
    # QP set up.
    m = length(mu_clin)
    Dm = diag(m+1) 
    Dm[nrow(Dm)-2, ncol(Dm)-2] <- Dm[nrow(Dm)-1, ncol(Dm)-1] <- 0.000001
    #Dm[nrow(Dm), ncol(Dm)] <- 1
    
    dv = rep(0,m+1)
    
    c_m1 <- c(1,1,0)
    c_m2 <- rbind(c(1,0,0),c(0,1,0))
    
    VY_diff <- (L_mu_mat[1,2] - L_mu_mat[-1,2])
    VZ_diff <- (L_mu_mat[1,3] - L_mu_mat[-1,3])
    V_diff <- cbind( VY_diff, VZ_diff, q)
    
    Am <- rbind(c_m1,c_m2, V_diff) 
    #Am <- rbind(c_m1, V_diff) 
    
    bv <- c(1,0.001,0.001,rep(0,k))
    #bv <- c(1,rep(0,k))
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
                                         M=Mk, lambda = lambda_est) }
    
    eta_opt_k <- show_condition(nmk(par = rep(0,5) , eta_func, control=list(maximize=TRUE))$par)
    if(length(eta_opt_k)==1){
      eta_opt_k <- show_condition(nmk(par = rep(0.01,5) , eta_func, control=list(maximize=TRUE))$par)
    }else if(length(eta_opt_k)==1){
      eta_opt_k <- show_condition(nmk(par = rep(0.02,5) , eta_func, control=list(maximize=TRUE))$par)
    }
    
    # Get decision rule.
    lin_modk <- as.vector(eta_opt_k %*% t(X))
    rulek <- trt_rule(lin_modk)
    
    VY_learner <- VY_est(rulek)
    VZ_learner <- VZ_tol_est(rulek, lambda_est)
    
    # Update data
    mu_learner <- c(VY_learner, -hinge(VZ_learner) )
    L_mu_mat <- rbind(L_mu_mat, c(-1,mu_learner) )
    
    
    #sum(rulek == A)
    # update iteration step.
    
    if(k >= 5 & abs(qk) <= eps ){break}
    
    k = k+1
    
  }
  
  M_est_vec[sim] <- Mk
  eta_mat_lambda[sim,] <-  eta_opt_k
  V_est_vec[sim] <-  mu_learner[1] - Mk*mu_learner[2]
  
  V_clin_vec[sim] <- mu_clin[1] - M_opt*mu_clin[2]
  
  VY_est_vec[sim] <- mu_learner[1]
  VZ_est_vec[sim] <- mu_learner[2]
  
  VY_clin_vec[sim] <- mu_clin[1]
  VZ_clin_vec[sim] <- mu_clin[2]
  
  PCD_vec[sim] <- sum(rulek==A_true)/n
  A_vec[sim] <- sum(A == A_true) / n
}

eta_bias <- apply(eta_mat_lambda,2,mean) - eta_opt

eta_sd <- apply(eta_mat_lambda,2,sd)

M_est_mean <- mean(M_est_vec)
M_est_sd <- sd(M_est_vec)

V_est_mean <- mean(V_est_vec)
V_est_sd <- sd(V_est_vec)

V_clin_mean <- mean(V_clin_vec)
V_clin_sd <- sd(V_clin_vec)

VY_est_mean <- mean(VY_est_vec)
VY_est_sd <- sd(VY_est_vec)

VY_clin_mean <- mean(VY_clin_vec)
VY_clin_sd <- sd(VY_clin_vec)

VZ_est_mean <- mean(VZ_est_vec)
VZ_est_sd <- sd(VZ_est_vec)


VZ_clin_mean <- mean(VZ_clin_vec)
VZ_clin_sd <- sd(VZ_clin_vec)

pcd_mean <- mean(PCD_vec)
pcd_sd <- sd(PCD_vec)

IRL_list <- list(
  eta_opt = eta_opt,
  eta_bias = eta_bias,
  eta_sd = eta_sd,
  M_est_mean = M_est_mean,
  M_est_sd = M_est_sd,
  V_est_mean = V_est_mean,
  V_est_sd = V_est_sd,
  V_clin_mean = V_clin_mean,
  V_clin_sd = V_clin_sd,
  VY_est_mean = VY_est_mean,
  VY_est_sd = VY_est_sd,
  VY_clin_mean = VY_clin_mean,
  VY_clin_sd = VY_clin_sd,
  VZ_est_mean = VZ_est_mean,
  VZ_est_sd = VZ_est_sd,
  VZ_clin_mean = VZ_clin_mean,
  VZ_clin_sd = VZ_clin_sd,
  pcd_mean = pcd_mean,
  pcd_sd = pcd_sd,
  PCD_vec = PCD_vec,
  M_est_vec = M_est_vec,
  A_mean = mean(A_vec),
  A_sd = sd(A_vec),
  clin_v_all=clin_v_all,
  eta_mat_est = eta_mat_lambda
  
)

capture.output(IRL_list, file = "IRL_known_lambda_settings4.txt")




