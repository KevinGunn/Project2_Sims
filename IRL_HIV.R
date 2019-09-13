setwd("~/EMR-research/Project 2 data")
HIVraw0<-read.csv("pmdata.csv",head=TRUE)
head(HIVraw0)
nrow(HIVraw0)

HIVraw<-HIVraw0
nrow(HIVraw)
table(HIVraw$drugc)
HIVnew<-subset(HIVraw,!is.na(drugc))
nrow(HIVnew)
HIVnew$Drug<-"PIs+NRTIs"
HIVnew$Drug[HIVnew$drugc==1]<-"NNRTIs+NRTIs"
table(HIVnew$Drug)
#####
library(rgenoud)
library(glmnet)
library(mstate)

#####
# Load data
mydata<-data.frame(time=HIVnew$FUtime,status=HIVnew$event,A=as.numeric(HIVnew$Drug=="NNRTIs+NRTIs"),CD4BL=HIVnew$CD4_BL,VLDBL=HIVnew$VLD_BL,ABNBL=HIVnew$ABNBL,
                    Age=scale(HIVnew$PatAge),Male=HIVnew$Male,Black=HIVnew$Black,GC=HIVnew$GC,Single=HIVnew$Single,Married=HIVnew$Married,Urban=HIVnew$urbann)

head(mydata)
#mydata<-subset(mydata0,!is.na(CD4BL))
n<-nrow(mydata)
summary(HIVnew$PatAge)
table(HIVnew$Male)
table(HIVnew$Black)
table(HIVnew$GC)


xx<-cbind(mydata$Age,mydata$Male,mydata$Black,mydata$GC,mydata$Single,mydata$Married,mydata$Urban)
colnames(xx)<-c("Age","Male","Black","GC","Single","Married","Urban")
fit<-cv.glmnet(x=xx,y=mydata$A,family="binomial",alpha=1)
coef(fit,s='lambda.min',exact=TRUE)


#############################################################################################
# function calculating the CIFs
calc.F<-function(beta,mydata){
  pix.hat<-predict(glm(A~ Age+Male+Black+GC+Single+Married+Urban,mydata,family=binomial("logit")),type="response")
  xb<-cbind(1,mydata$Age,mydata$Black,mydata$GC,mydata$Single,mydata$Married,mydata$Urban,mydata$Male)%*%beta
  
  #w.numer<-as.numeric(xb>=0)*mydata$A+as.numeric(xb<0)*(1-mydata$A)
  # smoothed version of weight
  n<-nrow(mydata)
  h<-sd(xb)*(n/4)^(-1/3)
  w.numer<-pnorm(xb/h)*mydata$A+(1-pnorm(xb/h))*(1-mydata$A)
  w.denom<-pix.hat*mydata$A+(1-pix.hat)*(1-mydata$A)
  wi<-w.numer/w.denom
  
  ord.t0<-sort(mydata$time)
  ord.d0<-as.numeric(mydata$status>0)[order(mydata$time)]
  ord.ebs<-mydata$status[order(mydata$time)]
  ord.wi<-wi[order(mydata$time)]
  
  tt<-ord.t0[ord.d0==1]
  SIt<-numeric()
  SIt[1]<-1
  for(k in 1:length(tt)){
    SIt[k+1]<-SIt[k]*(1-sum(ord.wi[ord.t0==tt[k]])/sum(ord.wi[ord.t0>=tt[k]]))
  }
  
  tt1<-ord.t0[ord.ebs==1]
  tt2<-ord.t0[ord.ebs==2]
  lamb1<-lamb2<-numeric()
  for(k in 1:length(tt1)) lamb1[k]<-sum(ord.wi[ord.t0==tt1[k]])/sum(ord.wi[ord.t0>=tt1[k]])
  for(k in 1:length(tt2)) lamb2[k]<-sum(ord.wi[ord.t0==tt2[k]])/sum(ord.wi[ord.t0>=tt2[k]])
  
  if(length(tt1)>0){ F1.tt0<-stepfun(tt1,cumsum(c(0,stepfun(tt,SIt)(tt1)*lamb1)))}else{F1.tt0<-0}
  if(length(tt2)>0){ F2.tt0<-stepfun(tt2,cumsum(c(0,stepfun(tt,SIt)(tt2)*lamb2)))}else{F2.tt0<-0}
  return(list(F1t0=F1.tt0,F2t0=F2.tt0,St0=stepfun(tt,SIt)))
}


# optim function
opt.beta1<-function(beta,mydata,tt0){
  cif.fit<-calc.F(beta,mydata)
  if(class(cif.fit$F1t0)[1]!="numeric"){F1.tt0<-cif.fit$F1t0(tt0)}else{F1.tt0<-0}
  if(class(cif.fit$F2t0)[1]!="numeric"){F2.tt0<-cif.fit$F2t0(tt0)}else{F2.tt0<-0}
  return(F1.tt0)
}

opt.beta2<-function(beta,mydata,tt0){
  cif.fit<-calc.F(beta,mydata)
  if(class(cif.fit$F1t0)[1]!="numeric"){F1.tt0<-cif.fit$F1t0(tt0)}else{F1.tt0<-0}
  if(class(cif.fit$F2t0)[1]!="numeric"){F2.tt0<-cif.fit$F2t0(tt0)}else{F2.tt0<-0}
  return(F2.tt0)
}

opt.beta<-function(beta,mydata,tt0=5,alp=0.1,M=1000){
  cif.fit<-calc.F(beta,mydata)
  if(class(cif.fit$F1t0)[1]!="numeric"){F1.tt0<-cif.fit$F1t0(tt0)}else{F1.tt0<-0}
  if(class(cif.fit$F2t0)[1]!="numeric"){F2.tt0<-cif.fit$F2t0(tt0)}else{F2.tt0<-0}
  return(F1.tt0+M*(F2.tt0-alp)*as.numeric(F2.tt0>alp))
}


# Condition handling code to avoid bad starting values in nmk.
show_condition <- function(code) {
  tryCatch(code,
           error = function(c) "error",
           warning = function(c) "warning",
           message = function(c) "message"
  )
}

# IRL QP program
QP_IRL <- function(k.num, tt0, eps, lambda, beta_initial){
  
  ##########################
  # Empty vectors/matrices to store values from each repetition.
  M_store <- rep(0,k.num)
  beta_store <- matrix(0,ncol=p,nrow=k.num)
  #########################
  
  # Clinicians Value function
  V_Y_clin <-  mean(Y_rc)
  V_Z_clin_tol <- - (mean(Z_ro) - lambda) 
  
  mu_clin <- c( V_Y_clin, -hinge(V_Z_clin_tol) )
  
  # Initialize decision rule
  lin_rule_coeff <- beta_initial     #c(0.01 , 0 , -0.5 , 0 , 0 , 0) 
  lin_mod1 <- as.vector(lin_rule_coeff %*% t(X.func))
  rule1 <- trt_rule(lin_mod1)
  
  #rule1
  print(sum(rule1 == HE_Cohort_Fluid_VP_sub$A))
  
  VY_learner <-  - VY_est(rule1, X.lm_int)
  VZ_learner <- VZ_tol_est(rule1, lambda, X.lm_int)
  
  mu_learner <- c(VY_learner, -hinge(VZ_learner) )
  
  
  # Create data matrix for QP program.
  mu_mat <- rbind(mu_clin, mu_learner)
  colnames(mu_mat) <- c("VY","VZ_tol")
  rownames(mu_mat) <- NULL
  
  Labels <- c(1,-1)
  L_mu_mat <- cbind(Labels,mu_mat)
  
  #set k to 1 to begin IRL QP.
  k=1
  
  # Dimensions of X
  p_beta <- dim(X.func)[2]
  
  while(k <= k.num){
    
    q = rep(1,k)
    
    # QP set up.
    m = length(mu_clin)
    Dm = diag(m+1) 
    Dm[nrow(Dm)-2, ncol(Dm)-2] <- Dm[nrow(Dm)-1, ncol(Dm)-1] <- 0.000001
    
    dv = rep(0,m+1)
    
    c_m1 <- c(1,1,0)
    c_m2 <- rbind(c(1,0,0),c(0,1,0))
    
    VY_diff <- (L_mu_mat[1,2] - L_mu_mat[-1,2])
    VZ_diff <- (L_mu_mat[1,3] - L_mu_mat[-1,3])
    V_diff <- cbind( VY_diff, VZ_diff, q)
    
    Am <- rbind(c_m1,c_m2, V_diff) 
    
    bv <- c(1,0.0001,0.0001,rep(0,k))
    
    # QP program
    sol <- solve.QP(Dm,dv,t(Am),bv,meq=1)
    
    Mk = sol$solution[2]/sol$solution[1]
    qk = sol$solution[3]
    
    # Store M values
    M_store[k] = Mk
    
    # Plot of Value functions.
    plot(L_mu_mat[,-1],col="blue", pch=19, main="Value Function Estimation",cex=1.2)
    points(x= L_mu_mat[1,2],y= L_mu_mat[1,3],bg="red",pch=22,cex=1.5) 
    # Hyperplane
    #abline( b=1/Mk , col="green", lty=4)
    
    
    # Obtain estimate of beta^opt.
  
    # Genetic Algorithm approach
    #tt0<-365*4
    for(j in 1:nrow(beta0)){
      fit<-try(optim(par=beta0[j,],fn=opt.beta,mydata=mydata,tt0=tt0,alp=alp,M=Mk),silent=TRUE)
      if(!is.character(fit)){
        if(fit$value<f3val){
          beta3<-fit$par/sqrt(sum(fit$par^2))
          beta3rec[l,]<-beta3
          f3val<-fit$value
          cat(paste("k=",k,"\n"))
          print(c(beta3,f3val))
        }
      }
    }
    fit<-try(genoud(opt.beta,nvars=npar,max=FALSE,starting.values=fit$par,max.generations=30,print.level=0,mydata=mydata,tt0=tt0,alp=alp,M=Mk),silent=TRUE)
    if(!is.character(fit)){
      if(fit$value<f3val){
        beta3<-fit$par/sqrt(sum(fit$par^2))
        beta3rec[l,]<-beta3
        f3val<-fit$value
        cat("rgenoud replace\n")
        print(c(beta3,f3val))
      }
    }
    
    # Store eta values
    eta_store[k,] <- eta_opt_k
    
    # Get decision rule.
    lin_modk <- as.vector(eta_opt_k %*% t(X.func))
    rulek <- trt_rule(lin_modk)
    
    VY_learner <- - VY_est(rulek,X.func)
    VZ_learner <- VZ_tol_est(rulek, lambda, X.func)
    
    # Update data
    mu_learner <- c(VY_learner, -hinge(VZ_learner) )
    L_mu_mat <- rbind(L_mu_mat, c(-1,mu_learner) )
    
    print(sum(rulek == HE_Cohort_Fluid_VP_sub$A) / dim(HE_Cohort_Fluid_VP_sub)[1])
    # update iteration step.
    
    if(k >= 5 & abs(qk) <= eps ){break}
    
    k = k+1
    
  }
  
  out_list <- list(M_store=M_store, 
                   eta_store = eta_store, 
                   eta_opt = eta_opt_k,
                   M_est = Mk,
                   eps_k = qk,
                   Value_mat = L_mu_mat,
                   trt_match = sum(rulek == HE_Cohort_Fluid_VP_sub$A) / dim(HE_Cohort_Fluid_VP_sub)[1],
                   decision_rule = rulek
  )
  
  return(out_list)
  
}

##################################################################################

# real data application
npar<-8
beta0<-matrix(0,2*npar,npar)
for(i in 1:npar){
  beta0[2*i-1,i]<-1
  beta0[2*i,i]<--1
}

mydata_r1 <- subset(mydata, status != 2)
CIF_1 <- Cuminc(mydata_r1$time, mydata_r1$status)

mydata_r2 <- subset(mydata, status != 1)
CIF_2 <- Cuminc(mydata_r2$time, mydata_r2$status)

CIFs_beta <- calc.F(beta0[1,], mydata)

#CIF for time=300
CIF_2[which.min(abs(300 - replace(CIF_2[,1], CIF_2[,1]>300, Inf))),3]
CIFs_beta$F2t0(300)

