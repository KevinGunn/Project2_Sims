setwd( "~/EMR-research/Project 2 data" )
sepdat<-read.csv( "rl_data_final_cont.csv" ,head=TRUE)
head(sepdat)
nrow(sepdat)

#####
library(rgenoud)
library(glmnet)
library(mstate)
library(quadprog)
library(cmprsk)
library(nloptr)
library(ggplot2)
library(gridExtra)
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
    # iv = 1 and vaso = 0
    iv_dose[i] = 1
  }
}

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

# Value Functions.

# Augmented IPWE
AIPWE_b <- function( Y, A, X, prop_score, eta){
  
  g_x_eta <- ifelse(eta %*% t(X) > 0, 1, 0)
  
  C_eta <- A*g_x_eta + (1-A)*( 1 - g_x_eta )
  
  prop_c <- prop_score*A + (1-prop_score)*( 1 - A )
  
  #dat <- data.frame(Y,A,X)
  #print(dat$X)
  mu <- glm(Y ~ X[,-1] + A*X[,-1], family = binomial(link = "logit"))
  
  A <- rep(1, length(Y))
  dat.m1 <- data.frame(X[,-1], A*X[,-1])
  A <- rep(0, length(Y))
  dat.m0 <- data.frame(X[,-1], A*X[,-1])
  
  mu_1 <- predict(mu, dat.m1)
  mu_0 <- predict(mu, dat.m0)
  
  m_x <- mu_1*g_x_eta + mu_0*( 1 - g_x_eta )
  
  AIPWE <- mean( ( as.vector( C_eta * Y ) / as.vector(prop_c) ) - 
                   ( as.vector( C_eta - as.vector(prop_c)) / as.vector(prop_c) )*as.vector(m_x) )
  
  return(AIPWE)
}

AIPWE_c <- function( Z, A, X, prop_score, eta){
  
  g_x_eta <- ifelse(eta %*% t(X) > 0, 1, 0)
  
  C_eta <- A*g_x_eta + (1-A)*( 1 - g_x_eta )
  
  prop_c <- prop_score*A + (1-prop_score)*( 1 - A )
  
  mu <- lm( Z ~ X[,-1] + A*X[,-1] )
  
  A <- rep(1, length(Z))
  dat.m1 <- data.frame(X[,-1], A*X[,-1])
  A <- rep(0, length(Z))
  dat.m0 <- data.frame(X[,-1], A*X[,-1])
  
  mu_1 <- predict(mu, dat.m1)
  mu_0 <- predict(mu, dat.m0)
  
  m_x <- mu_1*g_x_eta + mu_0*( 1 - g_x_eta )
  
  AIPWE <- mean( ( as.vector( C_eta * Z ) / as.vector(prop_c) ) - 
                   ( as.vector( C_eta - as.vector(prop_c)) / as.vector(prop_c) )*as.vector(m_x) )
  
  return(AIPWE)
}

# IPWE
IPWE_b <- function( Y, A, X, prop_score, eta){
  
  g_x_eta <- ifelse(eta %*% t(X) > 0, 1, 0)
  
  C_eta <- A*g_x_eta + (1-A)*( 1 - g_x_eta )
  
  prop_c <- prop_score*A + (1-prop_score)*( 1 - A )
  
  IPWE <- mean( as.vector( C_eta * Y ) / as.vector(prop_c ) )
  
  return(IPWE)
}

IPWE_c <- function( Z, A, X, prop_score, eta){
  
  g_x_eta <- ifelse(eta %*% t(X) > 0, 1, 0)
  
  C_eta <- A*g_x_eta + (1-A)*( 1 - g_x_eta )
  
  prop_c <- prop_score*A + (1-prop_score)*( 1 - A )
  
  IPWE <- mean( as.vector( C_eta * Z ) / as.vector(prop_c) )
  
  return(IPWE)
}

# NEED TO CHANGE IPWE IN OPT.ETA function.
opt.eta <- function(eta, Y, Z, A, X, prop_score, M, lambda){
  
  VY <- IPWE_b(Y, A, X, prop_score, eta)
  VZ <- IPWE_c(Z, A, X, prop_score, eta)
  
  return( -(VZ + M*(VY-lambda)*as.numeric(VY>lambda)) )
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
  V_Z_clin <-  as.numeric( mean(Z) )
  V_Y_clin_tol <- as.numeric( mean(Y) ) - lambda 
  
  mu_clin <- c( V_Z_clin, hinge(V_Y_clin_tol) )
  
  # Initialize decision rule    #c(0.01 , 0 , -0.5 , 0 , 0 , 0) 
  lin_mod1 <- as.vector(eta0 %*% t(X))
  rule1 <- trt_rule(lin_mod1)
  
  #rule1
  print( c("eta 0:", eta0) )
  print( c("rule 0", mean(rule1 == A)) )
  
  # Initial estimated rule
  VY_learner <- AIPWE_b(Y,A,X,prop_score,eta0)
  VZ_learner <- AIPWE_c(Z,A,X,prop_score,eta0)
  #VY_learner <- IPWE_b(Y,A,X,prop_score,eta0)
  #VZ_learner <- IPWE_c(Z,A,X,prop_score,eta0)
  
  mu_learner <- c( VZ_learner, hinge( VY_learner - lambda ) )
  
  # Create matrix for QP program.
  mu_mat <- rbind(mu_clin, mu_learner)
  colnames(mu_mat) <- c("VZ","VY_tol")
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
    Dm[nrow(Dm)-2, ncol(Dm)-2] <- Dm[nrow(Dm)-1, ncol(Dm)-1] <- 0.00000001
    
    dv = rep(0,m+1)
    
    c_m1 <- c(1,1,0)
    c_m2 <- rbind(c(1,0,0),c(0,1,0))
    
    VZ_diff <- (L_mu_mat[-1,2] - L_mu_mat[1,2] )
    VY_diff <- (L_mu_mat[-1,3] - L_mu_mat[1,3] )
    V_diff <- cbind( VZ_diff, VY_diff, q_coef)
    
    Am <- rbind(c_m1,c_m2, V_diff) 
    #Am <- rbind(c_m1,c_m2, V_diff[k,]) 
    
    bv <- c(1, 0.0000001, 0.0000001, rep(0,k))
    #bv <- c(1, 0.0000001, 0.0000001, rep(0,1))
    print(Am)
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
         xlab="Reward", ylab = "Proportion of Deaths in Hospital", cex=1.2, ylim = c(0,1), xlim = c(0,1))
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
    
    fval<-Inf
    for(j in 1:nrow(eta_start)){
      fit<-try(optim(par=eta_start[j,],fn=opt.eta,Y=Y,Z=Z,A=A,X=X,prop_score=prop_score,
                     lambda=lambda,M=Mk),silent=TRUE)
      if(!is.character(fit)){
        if(fit$value<fval){
          eta<- fit$par/sqrt(sum(fit$par^2))
          #etarec[l,]<-eta
          fval<-fit$value
          cat(paste("j=",j,"\n"))
          #print(c(eta,fval))
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
    VY_learner <- AIPWE_b(Y,A,X,prop_score,eta_k)
    VZ_learner <- AIPWE_c(Z,A,X,prop_score,eta_k)
    #VY_learner <- IPWE_b(Y,A,X,prop_score,eta_k)
    #VZ_learner <- IPWE_c(Z,A,X,prop_score,eta_k)
    
    # Update Value functions.
    mu_learner <- c(VZ_learner, hinge(VY_learner - lambda) )
    L_mu_mat <- rbind(L_mu_mat, c(-1,mu_learner) )
    
    #print(sum(rulek == data$A) / dim(data)[1])
    # update iteration step.
    
    if( (k >= 5 & abs(qk) <= eps) | (k >= 5 & all(eta_store[k-1,] == eta_k) ) ){break}
    
    k = k+1
    
  }
  
  out_list <- list(M_store=M_store, 
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

# Change rewards for those who survived/ did not survive in first bloc.
reward <- sepdat_iv_dose_b1$reward
reward[which(sepdat_iv_dose_b1$reward== 15)] <- 0
reward[which(sepdat_iv_dose_b1$reward== -15)] <- 0

# Recode died in hospital variable.
recode_dh <- sepdat_iv_dose_b1$died_in_hosp
recode_dh[which(sepdat_iv_dose_b1$died_in_hosp == 1)] <- 0
recode_dh[which(sepdat_iv_dose_b1$died_in_hosp == 0)] <- 1

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

# reduce covariates down to 10 if we set lambda to 0.03.# lambda = 0.005
lasso_prop <- glmnet(x=xvars,y=sepdat_iv_dose_b1$iv_dose, family="binomial",alpha=1, lambda=lambda.in)
coef(lasso_prop)
prop_fit <- predict(lasso_prop, type = "response", newx = xvars)
plot(lasso_prop, xvar = "lambda")

# covariates for final model
xvars_final <- xvars[, names(which(coef(lasso_prop)[-1,1] != 0)) ]
xvars_final <- cbind(1, xvars_final)

# real data application
npar<-dim(xvars_final)[2]
eta0<-matrix(0,2*npar,npar)
for(i in 1:npar){
  eta0[2*i-1,i]<- 1
  eta0[2*i,i]<- -1
}

#eta0_guess <- eta0[5,]
#eta0_guess <- eta0[4,]
#eta0_guess <- eta0[2,]
eta0_guess <- eta0[3,]

AIPWE_b(Y=sepdat_iv_dose_b1$died_in_hosp, A=sepdat_iv_dose_b1$iv_dose, X=xvars_final,prop_score = prop_fit,eta0_guess)
AIPWE_c(Z=sepdat_iv_dose_b1$reward, A=sepdat_iv_dose_b1$iv_dose, X=xvars_final,prop_score = prop_fit,eta0_guess)

IPWE_b(Y=sepdat_iv_dose_b1$died_in_hosp, A=sepdat_iv_dose_b1$iv_dose, X=xvars_final,prop_score = prop_fit,eta0_guess)
IPWE_c(Z=sepdat_iv_dose_b1$reward, A=sepdat_iv_dose_b1$iv_dose, X=xvars_final,prop_score = prop_fit,eta0_guess)

mean(sepdat_iv_dose_b1$died_in_hosp)
mean(sepdat_iv_dose_b1$reward)


# QP-IRL application.
AL_HIV <- QP_IRL( Y= recode_dh, Z=sepdat_iv_dose_b1$reward , A=sepdat_iv_dose_b1$iv_dose, 
                  X = xvars_final , prop_score = prop_fit , k.num=50, eps = 0.00001, 
                  lambda=0, eta0 = eta0_guess )


AIPWE_b(Y=sepdat_iv_dose_b1$died_in_hosp, A=sepdat_iv_dose_b1$iv_dose, X=xvars_final,prop_score = prop_fit,AL_HIV$eta_opt)
AIPWE_c(Z=sepdat_iv_dose_b1$reward, A=sepdat_iv_dose_b1$iv_dose, X=xvars_final,prop_score = prop_fit,AL_HIV$eta_opt)

IPWE_b(Y=sepdat_iv_dose_b1$died_in_hosp, A=sepdat_iv_dose_b1$iv_dose, X=xvars_final,prop_score = prop_fit,AL_HIV$eta_opt)
IPWE_c(Z=sepdat_iv_dose_b1$reward, A=sepdat_iv_dose_b1$iv_dose, X=xvars_final,prop_score = prop_fit,AL_HIV$eta_opt)

