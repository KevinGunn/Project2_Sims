###################################################
## Kevin Gunn                                    ##
## 12/4/2019                                     ##
## Project 2 - IRL                               ##
##                                               ##
###################################################

#####
library(rgenoud)
library(glmnet)
#library(mstate)
library(quadprog)
library(cmprsk)
library(nloptr)
library(ggplot2)
library(gridExtra)
library(MASS)
library(dfoptim)
#library(randomForest)
#####

###################################################
# Functions                                       #
###################################################
# Q_function for Value function
Q_coefs <- function(R, A, X){
  
  Qmod <- glm(R ~ X + A*X, family = binomial)
  
  coefs <- Qmod$coefficients
  
  return( coefs )
}


Q_V <- function(X, coefs, eta){
  
  g_x_eta <- ifelse(eta %*% t(X) > 0, 1, 0)
  
  xb <- as.vector(eta %*% t(X))
  g_x_eta <- ifelse(xb > 0, 1, 0)
  #print(dim(as.matrix(g_x_eta)*X))
  np <- dim(X)[2]
  
  X_FULL <- cbind(X, g_x_eta*X)
  
  V <- mean( round( 1 / ( 1 + exp(-( coefs %*% t(X_FULL) )) ) ) )
  
  return( V )
  
}

Q_opt.eta <- function(X,coef_Y, coef_Z, eta, lambda, M){
  
  g_x_eta <- ifelse(eta %*% t(X) > 0, 1, 0)
  
  xb <- as.vector(eta %*% t(X))
  g_x_eta <- ifelse(xb > 0, 1, 0)
  
  np <- dim(X)[2]
  
  X_FULL <- cbind(X, g_x_eta*X)
  
  VY <- mean( round( 1 / ( 1 + exp(-( coef_Y %*% t(X_FULL) )) ) ) )
  VZ <- mean( round( 1 / ( 1 + exp(-( coef_Z %*% t(X_FULL) )) ) ) )
  
  V_pen_est <- VY + M*hinge(VZ-lambda)
  
  return( V_pen_est )
  
}

# IPWE
IPWE <- function( R, A, X, prop_score, eta){
  
  #g_x_eta <- ifelse(eta %*% t(X) > 0, 1, 0)
  
  xb <- as.vector(eta %*% t(X))
  g_x_eta <- ifelse(xb > 0, 1, 0)
  
  C_eta <- A*g_x_eta + (1-A)*( 1 - g_x_eta )
  n<-length(A)
  #h<-sd(xb)*(n/4)^(-1/3)
  #C_eta <- pnorm(xb/h)*A + (1-pnorm(xb/h))*(1-A)
  
  #prop_c <- prop_score*A + (1-prop_score)*( 1 - A )
  prop_c <- prop_score*g_x_eta + (1-prop_score)*( 1 - g_x_eta )
  
  IPWE <- mean( as.vector( C_eta * R ) / as.vector(prop_c ) )
  
  return(IPWE)
}

# objective function.
opt.eta <- function(eta, Y, Z, A, X, prop_score, M, lambda){
  
  #VY <- AIPWE_b(Y, A, X, prop_score, eta)
  #VZ <- AIPWE_c(Z, A, X, prop_score, eta)
  
  VY <- IPWE(R=Y, A, X, prop_score, eta)
  VZ <- IPWE(R=Z, A, X, prop_score, eta)
  
  V_pen_est <- VY + M*hinge(VZ-lambda)
  
  return( V_pen_est )
}

# Obtain eta^{opt} and assume it is the clinicians decision rule.
clin_eta_opt <- function( eta, X, params_list, M, lambda){
  
  etax <- as.vector(eta %*% t(X))
  
  EY <-  mean( round( plogis( as.vector( params_list[[1]] %*% t(X) ) + 
                                trt_rule(etax) * as.vector( params_list[[2]] %*% t(X) ) ) ) )
  
  EZ_lambda <-  mean( round( plogis( as.vector( params_list[[3]] %*% t(X) ) + 
                                       trt_rule(etax) * as.vector( params_list[[4]] %*% t(X) ) ) ) ) - lambda
  
  # Penalized Value function.
  V_pen <- EY + M * hinge(EZ_lambda)
  
  #print(EY);print(hinge(EZ_lambda));print(V_pen)
  
  return(V_pen)
}

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

# QP IRL algorithm.
QP_IRL <- function(Y, Z, A, X, prop_score, k.num, eps, lambda, eta0){
  
  ##########################
  # Empty vectors/matrices to store values from each repetition.
  M_store <- rep(0,k.num)
  p = length(eta0)
  eta_store <- matrix(0,ncol=p,nrow=k.num)
  #########################
  
  # Clinicians Value function
  V_Y_clin <-  as.numeric( mean(Y) )
  V_Z_clin_tol <- as.numeric( mean(Z) ) - lambda 
  
  mu_clin <- c( V_Y_clin, hinge(V_Z_clin_tol) )
  
  # Initialize decision rule    #c(0.01 , 0 , -0.5 , 0 , 0 , 0) 
  lin_mod1 <- as.vector(eta0 %*% t(X))
  rule1 <- trt_rule(lin_mod1)
  
  #rule1
  print( c("eta 0:", eta0) )
  print( c("rule 0", mean(rule1 == A)) )
  
  # Initial estimated rule
  #VY_learner <- AIPWE_b(Y,A,X,prop_score,eta0)
  #VZ_learner <- AIPWE_c(Z,A,X,prop_score,eta0)
  #VY_learner <- IPWE(Y,A,X,prop_score,eta0)
  #VZ_learner <- IPWE(Z,A,X,prop_score,eta0)
  
  QY_coefs <- Q_coefs(Y, A, X[,-1])
  QZ_coefs <- Q_coefs(Z, A, X[,-1])
  
  VY_learner <- Q_V( X, QY_coefs, eta0 )
  VZ_learner <- Q_V( X, QZ_coefs, eta0 )
  
  mu_learner <- c( VY_learner, hinge( VZ_learner - lambda ) )
  
  # Create matrix for QP program.
  mu_mat <- rbind(mu_clin, mu_learner)
  colnames(mu_mat) <- c("VY","VZ_tol")
  rownames(mu_mat) <- NULL
  
  Labels <- c(1,-1)
  L_mu_mat <- cbind(Labels,mu_mat)
  
  #set k to 1 to begin IRL QP.
  k=1
  
  # Dimensions of X
  p_eta <- dim(X)[2]
  
  while(k <= k.num){
    
    q_coef = rep(1,k)
    
    # QP set up.
    m = length(mu_clin)
    Dm = diag(m+1) 
    Dm[nrow(Dm)-2, ncol(Dm)-2] <- Dm[nrow(Dm)-1, ncol(Dm)-1] <- 1e-16 #0.000001
    
    dv = rep(0,m+1)
    
    c_m1 <- c(1,1,0)
    c_m2 <- rbind(c(1,0,0),c(0,1,0))
    
    VY_diff <- (L_mu_mat[-1,2] - L_mu_mat[1,2] )
    VZ_diff <- (L_mu_mat[-1,3] - L_mu_mat[1,3] )
    V_diff <- cbind( VY_diff, VZ_diff, q_coef)
    
    #VY_diff <- (L_mu_mat[1,2] - L_mu_mat[-1,2] )
    #VZ_diff <- (L_mu_mat[1,3] - L_mu_mat[-1,3] )
    #V_diff <- cbind( VY_diff, VZ_diff, q_coef)
    
    Am <- rbind(c_m1,c_m2, V_diff) 
    
    bv <- c(1, 0.000001, 0.000001, rep(0,k))
    #bv <- c(1, 1e-16, 1e-16, rep(0,k))
    #bv <- c(1, 0, 0, rep(0,k))
    
    # QP program
    sol <- solve.QP(Dm,dv,t(Am),bv,meq=1)
    #print(sol)
    Mk = sol$solution[2]/sol$solution[1]
    qk = sol$solution[3]
    
    print(c(k,"Mk:",Mk)); print(c(k,"qk:",qk))
    print(c("weights", sol$solution[1],sol$solution[2]))
    
    # Store M values
    M_store[k] = Mk
    
    # Plot of Value functions.
    plot(L_mu_mat[-1,-1],col="blue", pch=19, main="Value Function Estimation",
         xlab="VY", ylab = "VZ", 
         cex=1.2, ylim = c(0,0.5), xlim = c(0,1))
    points(x= L_mu_mat[1,2],y= L_mu_mat[1,3],bg="red",pch=22,cex=1.5) 
    text(L_mu_mat[-1,-1], labels = seq(from=1, to = k) , pos = 4 )
    text(x = L_mu_mat[1,2], y= L_mu_mat[1,3], labels = "Clinicians" , pos = 1 )
    # Hyperplane
    #abline( b=1/Mk , col="green", lty=4)
    
    # Obtain estimate of eta^opt.
    
    # Genetic Algorithm approach
    #tt0<-365*4
    
    # real data application
    npar<-p
    eta_start<-matrix(0,2*npar,npar)
    for(i in 1:npar){
      eta_start[2*i-1,i]<- 1
      eta_start[2*i,i]<- -1
    }
    
    # record results in these matrices in case anything breaks.
    #numrec<-matrix(nrow=length(alps),ncol=4)
    #etarec<-matrix(nrow=length(alps),ncol=npar)
    #fvals<-matrix(nrow=length(alps),ncol=3)
    
    eta <- rep(0, npar)
    fval<-Inf
    for(j in 1:nrow(eta_start)){
      fit<-try(nmk(par=eta_start[j,],fn=Q_opt.eta, control = list(maximize=FALSE),
                   X=X,coef_Y=QY_coefs, coef_Z = QZ_coefs,
                   lambda=lambda,M=Mk),silent=TRUE)
      if(!is.character(fit)){
        if(fit$value<fval){
          eta <- fit$par/sqrt(sum(fit$par^2))
          #etarec[l,]<-eta
          fval<-fit$value
          cat(paste("j=",j,"\n"))
          print(c(eta,fval))
        }
      }
    }
    
    #fit<-try(genoud(opt.eta,nvars=npar,max=FALSE,starting.values=fit$par,max.generations=30,print.level=0,mydata=data,tt0=tt0.in,alp=lambda,M=Mk),silent=TRUE)
    #if(!is.character(fit)){
    #  if(fit$value<fval){
    #    eta<-fit$par/sqrt(sum(fit$par^2))
    #etarec[l,]<-eta
    #    fval<-fit$value
    #    cat("rgenoud replace\n")
    #    print(c(eta,fval))
    #  }
    #}
    
    #print(eta)
    eta_k <- eta
    print(c("eta k:", eta_k) )
    # Store eta values
    eta_store[k,] <- eta_k
    
    # Get decision rule.
    lin_modk <- as.vector(eta_k %*% t(X))
    #print(lin_modk)
    rulek <- trt_rule(lin_modk)
    
    print( c("rule_k",mean(rulek == A)) )
    
    # k estimated rule and value functions.
    #VY_learner <- AIPWE_b(Y,A,X,prop_score,eta_k)
    #VZ_learner <- AIPWE_c(Z,A,X,prop_score,eta_k)
    #VY_learner <- IPWE(Y,A,X,prop_score,eta_k)
    #VZ_learner <- IPWE(Z,A,X,prop_score,eta_k)
    VY_learner <- Q_V( X, QY_coefs, eta_k )
    VZ_learner <- Q_V( X, QZ_coefs, eta_k )
    
    # Update Value functions.
    mu_learner <- c(VY_learner, hinge(VZ_learner - lambda) )
    L_mu_mat <- rbind(L_mu_mat, c(-1,mu_learner) )
    
    #print(sum(rulek == data$A) / dim(data)[1])
    # update iteration step.
    
    if( (k >= 5 & abs(qk) <= eps) | (k >= 5 & all(eta_store[k-1,] == eta_k) ) ){break}
    
    k = k+1
    
  }
  
  out_list <- list(eta_initial = eta0,
                   M_store=M_store, 
                   eta_store = eta_store, 
                   eta_opt = eta_k,
                   M_est = Mk,
                   eps_k = qk,
                   Value_mat = L_mu_mat,
                   trt_match = sum(rulek == A) / length(A),
                   decision_rule = rulek
  )
  
  return(out_list)
  
}


#####################################################


#####################################################
# Provide Settings                                  #
#####################################################
set.seed(123)

p <- 5
n <- 10000
N <- 100000

#params

theta_0Y <- c(0,1,0,0,0,1)
#theta_0Y <- theta_0Y/sqrt(sum(theta_0Y^2))
theta_1Y <- c(0,1, 0.5, 0.5, 0.25, 1)
#theta_1Y <- theta_1Y/sqrt(sum(theta_1Y^2))
theta_0Z <- c(1,0,0,0,0,1)
#theta_0Z <- theta_0Z / sqrt(sum(theta_0Z^2))
theta_1Z <- c(-1,-0.75,0,-2,0,0)
#theta_1Z <- theta_1Z / sqrt(sum(theta_1Z^2))

params <- list(theta_0Y, theta_1Y, theta_0Z, theta_1Z)

# Generate Covariates 
X <- cbind(1, mvrnorm(n=N, mu=rep(0,p), diag(p) ) )

#####################################################################################


#####################################################################################
# Plot all values of M to get a reasonable value for their optimal 
# values.
#####################################################################################
M <- seq(0.1,1,0.1)
ml <- length(M)
npar<-p+1
lambda <- seq(0,0.4,0.2)

V <- rep(0,ml)
eta_opt_mat <- matrix(0,ncol=npar,nrow=ml)
EY_EZ <- matrix(0,ncol=2,nrow=ml)

eta_lambda_list <- list()
EY_EZ_lambda <- list()

for(l in lambda){
  
  for(i in 1:ml){
    
    # real data application
    eta_func <- function(eta){ clin_eta_opt(eta=eta, X=X, params_list = params,
                                            M=M[i], lambda = l)}
    
    eta_start<-matrix(0,2*npar,npar)
    for(k in 1:npar){
      eta_start[2*k-1,k]<- 1
      eta_start[2*k,k]<- -1
    }
    
    eta_opt <- rep(0,npar)
    
    fval<-Inf
    for(j in 1:nrow(eta_start)){
      fit<-try(nmk(par=eta_start[1,],fn=clin_eta_opt, control = list(maximize=FALSE),
                   X=X, params_list = params,M=M[i],lambda=l), silent = TRUE)
      if(!is.character(fit)){
        if(fit$value<fval){
          eta_opt <- fit$par/sqrt(sum(fit$par^2))
          #etarec[l,]<-eta
          fval<-fit$value
          cat(paste("j=",j,"\n"))
          print(c(eta_opt,fval))
        }
      }
    }
    
    
    print(eta_opt);print(eta_func(eta_opt))
    eta_opt_mat[i,] <- eta_opt
    V[i] <- eta_func(eta_opt)
    
    # Generate value functions for both outcomes
    EY_EZ[i,1] <- mean( round( plogis( as.vector( params[[1]] %*% t(X) ) + 
                                         trt_rule(as.vector(eta_opt %*% t(X))) * as.vector( params[[2]] %*% t(X) ) ) ) ) 
    
    
    EY_EZ[i,2] <-  hinge( mean( round( plogis( as.vector( params[[3]] %*% t(X) ) + 
                                                 trt_rule(as.vector(eta_opt %*% t(X))) * 
                                                 as.vector( params[[4]] %*% t(X) ) ) ) ) - l)
    
    
    cat("Value Functions:", EY_EZ[i,])
    
    
  }

  
  grid_vals <- cbind(M,V); colnames(grid_vals) <- c("M","V")
  gv_df <- as.data.frame(grid_vals)
  
  plot(x = gv_df[,1], y=gv_df[,2], main = "Value Function over M", col="blue")
  
  # list appending
  eta_lambda_list[[paste0("M", i, "lambda",l)]] <- eta_opt_mat
  EY_EZ_lambda[[paste0("M", i, "lambda",l)]] <- EY_EZ

}

# 
EY_EZ
eta_opt_mat
gv_df
M_opt = gv_df[3,1]; eta_opt <- eta_opt_mat[3,]


############################################################################################################
# number of simulations is set for 500.
set.seed(123)

num_sims <- 500
lambda_est <- 0

#Store output in following matrices/vectors.
eta_mat_lambda <- matrix(0, ncol=npar, nrow = num_sims)

M_est_vec <- rep(0, num_sims )
V_est_vec <- rep(0, num_sims )
V_clin_vec <- rep(0, num_sims )

VY_est_vec <- rep(0, num_sims )
VZ_est_vec <- rep(0, num_sims )


VY_clin_vec <- rep(0, num_sims )
VZ_clin_vec <- rep(0, num_sims )

clin_v_all <- matrix(0,ncol=2,nrow=num_sims)

clin_opt_match_vec <- rep(0, num_sims)

PCD_vec <- rep(0, num_sims )

AL_clin_vec <- rep(0, num_sims )

eta0 <- rep(0,npar)

for(sim in 1:num_sims){
  
  #### Generate all the data.
  # Generate Covariates Data 
  X <- cbind( 1, mvrnorm(n=n, mu=rep(0,p), diag(p)) )
  
  # 100% Correct decisions
  opt_model <- as.vector(eta_opt %*% t(X))
  true.prop <- plogis(opt_model)
  clin_model <- opt_model
  
  # Treatment
  A_opt <- trt_rule(clin_model)
  #runis <- runif(n,0,1)
  #A <- ifelse(runis < true.prop,1,0)
  A <- A_opt
  
  #prop <-cv.glmnet(x=X,y=A,family="binomial",alpha=1)
  #lasso_prop <- glmnet(x=X,y=A, family="binomial",alpha=1, lambda=prop$lambda.min)
  #prop_fit <- predict(lasso_prop, type = "response", newx = X)
  #hist(prop_fit)
  #A <- rbinom(n,1,prop_fit)
  
  #opt_trt <- trt_rule(opt_model)
  clin_opt_match <- mean(A == A_opt)
  clin_opt_match_vec[sim] <- clin_opt_match
  cat("A and opt Match",clin_opt_match, "\n")
  
  # Main Outcome of Interest
  #Y <- rbinom(n, 1, plogis( as.vector( theta_0Y %*% t(X) ) + A * as.vector( theta_1Y %*% t(X) ) ) )
  Y <- round( plogis( as.vector( theta_0Y %*% t(X) ) + A * as.vector( theta_1Y %*% t(X) ) ) )
  
  # Risk Outcome
  #Z <- rbinom(n, 1,  plogis( as.vector( theta_0Z %*% t(X) ) + A * as.vector( theta_1Z %*% t(X) ) ) )
  Z <- round( plogis( as.vector( theta_0Z %*% t(X) ) + A * as.vector( theta_1Z %*% t(X) ) ) )
  
  ####
  
  ####
  # Clinicians Value function
  V_Y_clin <-  mean(Y)
  V_Z_clin_tol <- mean(Z) - lambda_est  
  
  mu_clin <- c( V_Y_clin, hinge(V_Z_clin_tol) )
  
  clin_v_all[sim,] <- mu_clin
  
  # propensity score
  prop_score <- true.prop #prop_fit
  
  # QP-IRL application.
  AL_sim <- QP_IRL( Y= Y, Z=Z , A=A, X = X , prop_score = prop_score , k.num=50, eps = 0.00001, 
                    lambda=lambda_est, eta0 = eta0 )
  
  #eta_start<-matrix(0,2*npar,npar)
  #for(k in 1:npar){
  #  eta_start[2*k-1,k]<- 1
  #  eta_start[2*k,k]<- -1
  #}
  
  #eta_opt_guess <- rep(0,npar)
  #fval<-Inf
  #for(j in 1:nrow(eta_start)){
  #  fit<-try(nmk(par=eta_start[j,],fn=opt.eta, control = list(maximize=FALSE),
  #               Y=Y,Z=Z,A=A,X=X,prop_score=prop_score,
  #               lambda=lambda,M=M_opt),silent=TRUE)
  #  if(!is.character(fit)){
  #    if(fit$value<fval){
  #      eta_opt_guess <- fit$par/sqrt(sum(fit$par^2))
        #etarec[l,]<-eta
  #      fval<-fit$value
  #      cat(paste("j=",j,"\n"))
  #      print(c(eta_opt_guess,fval))
  #    }
  #  }
  #}
  
  
  M_est_vec[sim] <- AL_sim$M_est
  eta_mat_lambda[sim,] <-  AL_sim$eta_opt
  
  AL_iters <- dim(AL_sim$Value_mat)[1]
  V_est_vec[sim] <-  AL_sim$Value_mat[AL_iters,2] + AL_sim$M_est*AL_sim$Value_mat[AL_iters,3]
  
  V_clin_vec[sim] <- AL_sim$Value_mat[1,2] + AL_sim$M_est*AL_sim$Value_mat[1,3]
  
  VY_est_vec[sim] <- AL_sim$Value_mat[AL_iters,2]
  VZ_est_vec[sim] <- AL_sim$Value_mat[AL_iters,3]
  
  VY_clin_vec[sim] <- AL_sim$Value_mat[1,2]
  VZ_clin_vec[sim] <- AL_sim$Value_mat[1,3]
  
  AL_clin_vec[sim] <- AL_sim$trt_match
  PCD_vec[sim] <- mean(AL_sim$decision_rule == A_opt)
  
}

eta_mean <- apply(eta_mat_lambda,2,mean)

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

AL_clin_match_mean <- mean(AL_clin_vec)
AL_clin_match_sd <- sd(AL_clin_vec)

PCD_mean <- mean(PCD_vec)
PCD_sd <- sd(PCD_vec)


clin_opt_match_mean <- mean(clin_opt_match_vec)
clin_opt_match_sd <- sd(clin_opt_match_vec)

IRL_list <- list(
  eta_opt = eta_opt,
  eta_mean = eta_mean,
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
  pcd_mean = PCD_mean,
  pcd_sd = PCD_sd,
  AL_clin_match_mean = AL_clin_match_mean,
  AL_clin_match_sd = AL_clin_match_sd,
  pcd_vec = PCD_vec,
  M_est_vec = M_est_vec,
  clin_v_all=clin_v_all,
  clin_opt_match_mean = clin_opt_match_mean, 
  clin_opt_match_sd = clin_opt_match_sd,
  eta_mat_est = eta_mat_lambda
  
)

capture.output(IRL_list, file = "IRL_min_func_smallM.txt")

# Things to do:
# record clinician  pcd
# record algorithm matches pcd
