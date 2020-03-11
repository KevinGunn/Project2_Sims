setwd( "~/EMR-research/Project 2 data" )
sepdat<-read.csv( "rl_data_final_cont.csv" ,head=TRUE)
head(sepdat)
nrow(sepdat)


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

# look at the available variables.
colnames(sepdat)

######################################################################################################

# Extract subjects who only received IV_fluids or no treatment.
sepdat_iv <- sepdat[ which(sepdat$vaso_input==0), ] 


######################################################################################################
# Plots of IV treatment trajectory.

thedat <- within(sepdat_iv, {  bloc <- factor(bloc)
iv_input <- factor(iv_input,levels=0:4)
id <- factor(icustayid)
})

#   First create the basic plot object
pp <- ggplot(thedat,aes(x=bloc,y=iv_input,group=id))

#Average Dosage Trajectory using all 20955 subjects.
pp + stat_summary(aes(group = 1,color= "red" ),
                  geom =  "line" , fun.y = mean, size = 1.5) +
  xlab( "Bloc" ) + ylab( "Dosage" ) + theme(legend.position =  "none" ) +
  ggtitle( "Average Trajectory for Subjects Receiving only IV Fluid" )
#dev.off()


# Random subset of patients
sepdat_5 <- sepdat_iv[ which( sepdat_iv$icustayid %in% sepdat_iv$icustayid[sample(nrow(sepdat_iv),5)] ) ,  ]

thedat <- within(sepdat_5, {  bloc <- factor(bloc)
iv_input <- factor(iv_input,levels=0:4)
id <- factor(icustayid)
})

#   First create the basic plot object
pp2 <- ggplot( thedat,aes(x=bloc,y=iv_input,group=id, color=id, shape=id ))

pp2 + geom_line() + geom_point(shape=18) + xlab( "Bloc" ) + ylab( "Dosage" ) +  
  theme(legend.position =  "none" )

##########################################################################################################

# Count number of subjects in high dose and low dose treatment groups at bloc 1.
## Pull out first bloc for each subject.

sepdat_iv_1 <- subset(sepdat_iv,bloc==1)

length(which(sepdat_iv_1$iv_input >= 3))

length(which(sepdat_iv_1$iv_input < 3))



# Count number of subjects that change treatment throughout stay.
# Ex. high dose to low dose or low dose to high dose.

# Step 1. - rewrite iv trt as low dose and high dose

iv_dose <- rep(0,dim(sepdat_iv)[1] )

for( i in 1:dim(sepdat_iv)[1] ){
  
  if( sepdat_iv$iv_input[i] < 3 ){ 
    # No trt
    iv_dose[i] = 0
    
  } else if ( sepdat_iv$iv_input[i] >= 3 ){
    #iv = 1 and vaso = 0
    iv_dose[i] = 1
  }
}



#iv_dose <- rep(0,dim(sepdat_iv)[1] )

#for( i in 1:dim(sepdat_iv)[1] ){

#  if( sepdat_iv$iv_input[i] == 0 ){ 
# No trt
#    iv_dose[i] = 0

#  } else if ( sepdat_iv$iv_input[i] > 0 ){
# iv = 1 and vaso = 0
#    iv_dose[i] = 1
#  }
#}

sepdat_iv_dose <- data.frame(sepdat_iv, iv_dose)

# Step 2. - For each block count number of subjects in each trt group by storing in matrix.

trt_bloc_mat <- matrix(0,ncol=3,nrow=20)
colnames(trt_bloc_mat) <- c( "bloc" ,  "Low Dose" ,  "High Dose" )

blocs <- sort(unique(sepdat_iv_dose$bloc))

for( i in blocs ){
  
  sepdat_iv_dose_bloc <- subset(sepdat_iv_dose, bloc==i)
  bt <- table(sepdat_iv_dose_bloc$iv_dose)
  trt_bloc_mat[i,] <- c( i, bt[1], bt[2] )
  
}


plot(trt_bloc_mat[,2], ylim = c(0,18000), xlab= "Blocks" , ylab= "Number of Patients" )
lines(trt_bloc_mat[,2])
points(trt_bloc_mat[,3], col=2)
lines(trt_bloc_mat[,3],col=2)
# Add a legend
legend( "topright" , legend=c( "Low Dose" ,  "High Dose" ),
        col=c( "black" ,  "red" ), lty=1, cex=0.8)

# Step 3. - Make a table.
library(xtable)

xtable(trt_bloc_mat)

#######################################################################################################

# Analyze data only using Block 1.
sepdat_iv_dose_b1 <- subset(sepdat_iv_dose,bloc==1)

# Propensity Score Analysis.
prop_fit <- glm(iv_dose ~ gender + age + elixhauser + Weight_kg + GCS + HR + SysBP + MeanBP + DiaBP + RR +                                         
                  SpO2 + Temp_C + FiO2_1 + Potassium + Sodium + Chloride + Glucose + BUN + Creatinine +                                
                  Magnesium + Calcium + Ionised_Ca + CO2_mEqL + SGOT + SGPT + Total_bili + Albumin + Hb +                                         
                  WBC_count + Platelets_count + PTT + PT + INR + Arterial_pH + paO2 + paCO2 + Arterial_BE +                                
                  Arterial_lactate + HCO3 + mechvent + Shock_Index + PaO2_FiO2 + SOFA + SIRS, 
                data = sepdat_iv_dose_b1, family=binomial(link = "logit")              
)

#summary(prop_fit)

names( which( coef(summary(prop_fit))[,4] < 0.05 ) ) 

names( which( coef(summary(prop_fit))[,4] < 0.000001 ) ) 

xtable(summary(prop_fit))


mean( sepdat_iv_dose_b1$died_in_hosp )

mean( sepdat_iv_dose_b1$reward )

#######################################################################################################

# Code to implement IRL algorithm.

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
  
  V <- mean( ( 1 / ( 1 + exp(-( coefs %*% t(X_FULL) )) ) ) )
  
  return( V )
  
}

Q_opt.eta <- function(X,coef_Y, coef_Z, eta, lambda, M){
  
  g_x_eta <- ifelse(eta %*% t(X) > 0, 1, 0)
  
  xb <- as.vector(eta %*% t(X))
  g_x_eta <- ifelse(xb > 0, 1, 0)
  
  np <- dim(X)[2]
  
  X_FULL <- cbind(X, g_x_eta*X)
  
  VY <- mean( ( 1 / ( 1 + exp(-( coef_Y %*% t(X_FULL) )) ) ) )
  VZ <- mean( ( 1 / ( 1 + exp(-( coef_Z %*% t(X_FULL) )) ) ) )
  
  V_pen_est <- VY + M*hinge(VZ-lambda)
  
  return( V_pen_est )
  
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
QP_IRL <- function(Y, Z, A, X, k.num, eps, lambda, eta0){
  
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
    #Dm[nrow(Dm), ncol(Dm)] <- 1e-16 #0.000001
    
    dv = rep(0,m+1)
    
    c_m1 <- c(1,1,0)
    c_m2 <- rbind(c(1,0,0),c(0,1,0))
    
    #VY_diff <- (L_mu_mat[-1,2] - L_mu_mat[1,2] )
    #VZ_diff <- (L_mu_mat[-1,3] - L_mu_mat[1,3] )
    #V_diff <- cbind( VY_diff, VZ_diff, q_coef)
    
    VY_diff <- (L_mu_mat[1,2] - L_mu_mat[-1,2] )
    VZ_diff <- (L_mu_mat[1,3] - L_mu_mat[-1,3] )
    V_diff <- cbind( VY_diff, VZ_diff, q_coef)
    
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
         cex=1.2, ylim = c(0.7,1), xlim = c(0.7,1))
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
      fit<-try(nmk(par=eta_start[j,],fn=Q_opt.eta, control = list(maximize=TRUE),
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
    VY_learner <- Q_V( X, QY_coefs, eta_k )
    VZ_learner <- Q_V( X, QZ_coefs, eta_k )
    
    # Update Value functions.
    mu_learner <- c(VY_learner, hinge(VZ_learner - lambda) )
    L_mu_mat <- rbind(L_mu_mat, c(-1,mu_learner) )
    
    #print(sum(rulek == data$A) / dim(data)[1])
    # update iteration step.
    
    if( (k >= 10 & abs(qk) <= eps) | (k >= 10 & all(eta_store[k-1,] == eta_k) ) ){break}
    
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


#########################################################################################################

# Real Data Application.

# Recode died in hospital variable.
recode_dh <- sepdat_iv_dose_b1$died_in_hosp
recode_dh[which(sepdat_iv_dose_b1$died_in_hosp == 1)] <- 0
recode_dh[which(sepdat_iv_dose_b1$died_in_hosp == 0)] <- 1

recode_d90 <- sepdat_iv_dose_b1$mortality_90d
recode_d90[which(sepdat_iv_dose_b1$mortality_90d == 1)] <- 0
recode_d90[which(sepdat_iv_dose_b1$mortality_90d == 0)] <- 1

recode_d48 <- sepdat_iv_dose_b1$died_within_48h_of_out_time
recode_d48[which(sepdat_iv_dose_b1$died_within_48h_of_out_time == 1)] <- 0
recode_d48[which(sepdat_iv_dose_b1$died_within_48h_of_out_time == 0)] <- 1

# Covariates for glmnet
xvars_noscale <- cbind(sepdat_iv_dose_b1$gender, sepdat_iv_dose_b1$age, sepdat_iv_dose_b1$elixhauser, 
                       sepdat_iv_dose_b1$Weight_kg, sepdat_iv_dose_b1$GCS, sepdat_iv_dose_b1$HR, sepdat_iv_dose_b1$SysBP,
                       sepdat_iv_dose_b1$MeanBP, sepdat_iv_dose_b1$DiaBP, sepdat_iv_dose_b1$RR,                                         
                       sepdat_iv_dose_b1$SpO2, sepdat_iv_dose_b1$Temp_C, sepdat_iv_dose_b1$FiO2_1, 
                       sepdat_iv_dose_b1$Potassium, sepdat_iv_dose_b1$Sodium, sepdat_iv_dose_b1$Chloride,
                       sepdat_iv_dose_b1$Glucose, sepdat_iv_dose_b1$BUN, sepdat_iv_dose_b1$Creatinine,                                
                       sepdat_iv_dose_b1$Magnesium, sepdat_iv_dose_b1$Calcium, sepdat_iv_dose_b1$Ionised_Ca,
                       sepdat_iv_dose_b1$CO2_mEqL, sepdat_iv_dose_b1$SGOT, sepdat_iv_dose_b1$SGPT, sepdat_iv_dose_b1$Total_bili, sepdat_iv_dose_b1$Albumin, 
                       sepdat_iv_dose_b1$Hb, sepdat_iv_dose_b1$WBC_count, sepdat_iv_dose_b1$Platelets_count,
                       sepdat_iv_dose_b1$PTT, sepdat_iv_dose_b1$PT, sepdat_iv_dose_b1$INR, sepdat_iv_dose_b1$Arterial_pH,
                       sepdat_iv_dose_b1$paO2, sepdat_iv_dose_b1$paCO2, sepdat_iv_dose_b1$Arterial_BE, sepdat_iv_dose_b1$Arterial_lactate,
                       sepdat_iv_dose_b1$HCO3, sepdat_iv_dose_b1$mechvent, sepdat_iv_dose_b1$Shock_Index,
                       sepdat_iv_dose_b1$PaO2_FiO2, sepdat_iv_dose_b1$SOFA, sepdat_iv_dose_b1$SIRS)

xvars <- cbind(xvars_noscale[, c(1,40)], log(xvars_noscale[, 43 ]+1), log( xvars_noscale[, 44] + 1),
               apply( xvars_noscale[, -c(1,40,43,44)], 2 , scale ) )

colnames(xvars) <- c( "gender", "mechvent", "SOFA", "SIRS", "age", "elixhauser", "Weight_kg", "GCS",
                      "HR", "SysBP", "MeanBP", "DiaBP", "RR", "SpO2", "Temp_C", "FiO2_1", "Potassium", 
                      "Sodium", "Chloride", "Glucose", "BUN", "Creatinine", "Magnesium", "Calcium",
                      "Ionised_Ca", "CO2_mEqL", "SGOT", "SGPT", "Total_bili", "Albumin", "Hb",                                         
                      "WBC_count", "Platelets_count", "PTT", "PT", "INR", "Arterial_pH", "paO2", "paCO2",
                      "Arterial_BE", "Arterial_lactate", "HCO3", "Shock_Index", "PaO2_FiO2")

# Propensity Score with glmnet.
fit <- cv.glmnet(x=xvars,y=sepdat_iv_dose_b1$iv_dose, family="binomial",alpha=1)
lambda.in <- fit$lambda.1se
# Covariates of interest.
coef(fit,s='lambda.1se',exact=TRUE)

# Covariates of interest.
# 
#lasso_prop <- glmnet(x=xvars,y=sepdat_iv_dose_b1$iv_dose, family="binomial",alpha=1, lambda=lambda.in)
#prop_fit <- predict(lasso_prop, type = "response", newx = xvars)

#lambda = 0.035 sets the number of variables to about 8.
# Use BIC to pick variables.
n_vars_out <- rep(0,10)
lv <- seq(0.01,0.05, by=0.005)
for(ll in 1:length(lv)){
  lasso_prop <- glmnet(x=xvars,y=sepdat_iv_dose_b1$iv_dose, family="binomial",alpha=1, lambda=lv[ll])
  prop_fit <- predict(lasso_prop, type = "response", newx = xvars)
  
  xvars_final <- xvars[, names(which(coef(lasso_prop)[-1,1] != 0)) ]
  xvars_final <- cbind(1, xvars_final)
  
  # real data application
  npar<-dim(xvars_final)[2]
  n_vars_out[ll] <- npar
}

print(n_vars_out)

# Reduce down to 6 variables.
lasso_prop <- glmnet(x=xvars,y=sepdat_iv_dose_b1$iv_dose, family="binomial",alpha=1, lambda=lv[5])
prop_fit <- predict(lasso_prop, type = "response", newx = xvars)

xvars_final <- xvars[, names(which(coef(lasso_prop)[-1,1] != 0)) ]
xvars_final <- cbind(1, xvars_final)

npar<-dim(xvars_final)[2]

eta0_1 <- rep(0,npar)
eta0_2 <- rep(-1,npar)
eta0_3 <- rep(1,npar)

# QP-IRL application.
AL_sep_90 <- QP_IRL( Y= recode_dh, Z=recode_d90 , A=sepdat_iv_dose_b1$iv_dose, 
                   X = xvars_final , k.num=50, eps = 0.00001, 
                   lambda=0, eta0 = eta0_1 )


#capture.output(AL_sep_90,file = "AL_sep_90_lambda_035.txt")

# QP-IRL application.
AL_sep_48 <- QP_IRL( Y= recode_dh, Z=recode_d48 , A=sepdat_iv_dose_b1$iv_dose, 
                  X = xvars_final , k.num=50, eps = 0.00001, 
                  lambda=0, eta0 = eta0_1 )

capture.output(AL_sep_48,file = "AL_sep_48_start1.txt")

##############################################################################

# OTD by optimizing Y and Z seperately.

QY_coefs <- Q_coefs(recode_dh, sepdat_iv_dose_b1$iv_dose,  xvars_final[,-1])
QZ_coefs <- Q_coefs(recode_d48, sepdat_iv_dose_b1$iv_dose,  xvars_final[,-1])

# QP IRL algorithm.

npar<-dim(xvars_final)
eta_start<-matrix(0,2*npar,npar)
for(i in 1:npar){
    eta_start[2*i-1,i]<- 1
    eta_start[2*i,i]<- -1
  }
    

fval<-Inf
for(j in 1:nrow(eta_start)){
      fit<-try(nmk(par=eta_start[j,],fn=Q_V, control = list(maximize=TRUE),
                   X= xvars_final,coefs=QY_coefs),silent=TRUE)
      if(!is.character(fit)){
        if(fit$value<fval){
          eta_opt_Y <- fit$par/sqrt(sum(fit$par^2))
          #etarec[l,]<-eta
          fval<-fit$value
          cat(paste("j=",j,"\n"))
          print(c(eta_opt_Y,fval))
        }
      }
    }
    
fval<-Inf
for(j in 1:nrow(eta_start)){
  fit<-try(nmk(par=eta_start[j,],fn=Q_V, control = list(maximize=TRUE),
               X= xvars_final,coefs=QZ_coefs),silent=TRUE)
  if(!is.character(fit)){
    if(fit$value<fval){
      eta_opt_Z <- fit$par/sqrt(sum(fit$par^2))
      #etarec[l,]<-eta
      fval<-fit$value
      cat(paste("j=",j,"\n"))
      print(c(eta_opt_Z,fval))
    }
  }
}

A_opt_Y <- trt_rule(eta_opt_Y %*% t(xvars_final))
A_opt_Z <- trt_rule(eta_opt_Z %*% t(xvars_final))

mean(A_opt_Y == A_opt_Z)
# Agree 72% of the time.
##############################################################################
# Random Forest approach
##############################################################################

#####
library(randomForestSRC)
#####

###################################################
# Functions                                       #
###################################################
# Q_function for Value function

Q_V <- function(pred_mat, X, eta){
  
  #Pred mat to vectors
  pred0 <- as.vector(pred_mat[,1])
  pred1 <- as.vector(pred_mat[,2])
  
  #Decision rule
  xb <- as.vector(eta %*% t(X))
  g_x_eta <- ifelse(xb > 0, 1, 0)
  
  #Value for R.
  V <- mean( (1-g_x_eta)*pred0 + g_x_eta*pred1 )
  
  return( V )
  
}

Q_preds <- function(R, X, A){
  
  # Stratifuy by trt when using random forests.
  X0 <- X[which(A==0), ]
  R0 <- as.factor(R[which(A==0)])
  d0 <- data.frame(R0 = R0, X = X0)
  
  X1 <- X[which(A==1), ]
  R1 <- as.factor(R[which(A==1)])
  d1 <- data.frame(R1 = R1, X = X1)
  
  control_synth <- synthetic(R0 ~. , data = d0 )
  trt_synth  <- synthetic(R1 ~. , data = d1 )
  
  ##get appropriate out of bag predictions to estimate Q1 and Q0
  ##oob control pred (control pred using control data)
  c_object <- control_synth$rfSyn
  c_oob <-  c_object$predicted.oob[,2]
  
  ##control pred using trt data
  c_tdata <- synthetic(object = control_synth, newdata = d1)$rfSynPred$predicted[,2]
  
  ##oob trt pred (trt pred using trt data)
  t_object <- trt_synth$rfSyn
  t_oob <- t_object$predicted.oob[,2]
  
  ##trt pred using control data
  t_cdata <- synthetic(object = trt_synth, newdata = d0)$rfSynPred$predicted[,2]
  
  #Combine matrix
  X_all <- as.matrix(rbind(d0[,-1],d1[,-1]))
  
  #Combine preds
  pred0 <- c(c_oob , c_tdata)
  pred1 <- c(t_cdata , t_oob)
  
  pred_mat <- cbind(pred0, pred1)
  colnames(pred_mat) <- c("A0","A1")
  
  out_list <- list(pred_mat = pred_mat, X= X_all)
  
  return(out_list)
}

Q_preds_LR <- function(R, X, A){
  
  # Stratifuy by trt when using random forests.
  X0 <- X[which(A==0), ]
  R0 <- as.factor(R[which(A==0)])
  d0 <- data.frame(R = R0, X = X0)
  
  X1 <- X[which(A==1), ]
  R1 <- as.factor(R[which(A==1)])
  d1 <- data.frame(R = R1, X = X1)
  
  #Combine matrix
  X_all <- rbind(d0[,-1],d1[,-1])
  colnames(d1) <- colnames(d0)
  
  #################################################################
  
  # RF for trt 0
  LR_R0 <- glm(R ~ ., data=d0[,-2], family=binomial(link = "logit"))
  pred0 <- predict(LR_R0, X_all, type="response")
  
  # RF for trt 1
  LR_R1 <- glm(R ~ ., data=d1[,-2], family=binomial(link = "logit"))
  pred1 <- predict(LR_R1, X_all, type="response")
  
  #Combine matrix
  X_all <- as.matrix(X_all)
  #Combine preds
  pred_mat <- cbind(pred0, pred1)
  colnames(pred_mat) <- c("A0","A1")
  
  out_list <- list(pred_mat = pred_mat, X = X_all)
  
  return( out_list )
  
}

Q_opt.eta <- function(pred_mat_Y,pred_mat_Z, X, eta, lambda, M){
  
  #Decision rule
  xb <- as.vector(eta %*% t(X))
  g_x_eta <- ifelse(xb > 0, 1, 0)
  
  #Pred mat to vectors
  predY0 <- as.vector(pred_mat_Y[,1])
  predY1 <- as.vector(pred_mat_Y[,2])
  
  predZ0 <- as.vector(pred_mat_Z[,1])
  predZ1 <- as.vector(pred_mat_Z[,2])
  
  #Value for Y and Z.
  VY <- mean( ( 1- g_x_eta )*predY0 + g_x_eta*predY1 )
  VZ <- mean( ( 1- g_x_eta )*predZ0 + g_x_eta*predZ1 )
  
  V_pen_est <- VY + M*hinge(VZ-lambda)
  
  return( V_pen_est )
  
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
QP_IRL <- function(Y, Z, A, X, k.num, eps, lambda, eta0){
  
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
  
  # Q function estimates for Y and Z.
  QY_list <- Q_preds(Y, X, A)
  QZ_list <- Q_preds(Z, X, A)
  
  QY_list_LR <- Q_preds_LR(Y, X, A)
  QZ_list_LR <- Q_preds_LR(Z, X, A)
  
  QY1_MSE <- mean( (QY_list$pred_mat[,2] - QY_list_LR$pred_mat[,2])^2  )
  QY0_MSE <- mean( (QY_list$pred_mat[,1] - QY_list_LR$pred_mat[,1])^2  )
  
  QZ1_MSE <- mean( (QZ_list$pred_mat[,2] - QZ_list_LR$pred_mat[,2])^2  )
  QZ0_MSE <- mean( (QZ_list$pred_mat[,1] - QZ_list_LR$pred_mat[,1])^2  )
  
  # Initial estimated rule
  VY_learner <- Q_V( QY_list$pred_mat, QY_list$X, eta0 )
  VZ_learner <- Q_V( QZ_list$pred_mat, QZ_list$X, eta0 )
  
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
    #Dm[nrow(Dm), ncol(Dm)] <- 1e-16 #0.000001
    
    dv = rep(0,m+1)
    
    c_m1 <- c(1,1,0)
    c_m2 <- rbind(c(1,0,0),c(0,1,0))
    
    #VY_diff <- (L_mu_mat[-1,2] - L_mu_mat[1,2] )
    #VZ_diff <- (L_mu_mat[-1,3] - L_mu_mat[1,3] )
    #V_diff <- cbind( VY_diff, VZ_diff, q_coef)
    
    VY_diff <- (L_mu_mat[1,2] - L_mu_mat[-1,2] )
    VZ_diff <- (L_mu_mat[1,3] - L_mu_mat[-1,3] )
    V_diff <- cbind( VY_diff, VZ_diff, q_coef)
    
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
      fit<-try(nmk(par=eta_start[j,],fn=Q_opt.eta, control = list(maximize=TRUE),
                   pred_mat_Y=QY_list$pred_mat, pred_mat_Z=QZ_list$pred_mat, X=QY_list$X,
                   lambda=lambda,M=Mk),silent=TRUE)
      if(!is.character(fit)){
        if(fit$value<=fval){
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
    VY_learner <- Q_V( QY_list$pred_mat, QY_list$X, eta_k )
    VZ_learner <- Q_V( QZ_list$pred_mat, QZ_list$X, eta_k )
    
    # Update Value functions.
    mu_learner <- c(VY_learner, hinge(VZ_learner - lambda) )
    L_mu_mat <- rbind(L_mu_mat, c(-1,mu_learner) )
    
    #print(sum(rulek == data$A) / dim(data)[1])
    # update iteration step.
    
    if( (k >= 10 & qk <= eps) | (k >= 10 & all(eta_store[k-1,] == eta_k) ) ){break}
    
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
                   decision_rule = rulek,
                   QY1_MSE = QY1_MSE,
                   QY0_MSE = QY0_MSE,
                   QZ1_MSE = QZ1_MSE,
                   QZ0_MSE = QZ0_MSE
  )
  
  return(out_list)
  
}

# QP-IRL application.
AL_sep_48_RF <- QP_IRL( Y= recode_dh, Z=recode_d48 , A=sepdat_iv_dose_b1$iv_dose, 
                     X = xvars_final , k.num=50, eps = 0.00001, 
                     lambda=0, eta0 = eta0_2 )

capture.output(AL_sep_48_RF,file = "AL_sep_48_RF.txt")

###############################################################################
