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
library(quadprog)

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
calc.F<-function(eta,mydata){
  pix.hat<-predict(glm(A~ Age+Male+Black+GC+Single+Married+Urban,mydata,family=binomial("logit")),type="response")
  xb<-cbind(1,mydata$Age,mydata$Black,mydata$GC,mydata$Single,mydata$Married,mydata$Urban,mydata$Male)%*%eta
  
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
opt.eta<-function(eta,mydata,tt0=5,alp=0.1,M=1000){
  cif.fit<-calc.F(eta,mydata)
  if(class(cif.fit$F1t0)[1]!="numeric"){F1.tt0<-cif.fit$F1t0(tt0)}else{F1.tt0<-0}
  if(class(cif.fit$F2t0)[1]!="numeric"){F2.tt0<-cif.fit$F2t0(tt0)}else{F2.tt0<-0}
  return( F1.tt0+M*(F2.tt0-alp)*as.numeric(F2.tt0>alp) )
}

# Penalty Function and linear decision rule functions:

norm <- function(x) sqrt(sum(x^2))

hinge <- function(s){ ifelse(s>0,s,0) } 

trt_rule <- function(z){ ifelse(z>0,1,0) }


# Clinician's CIF.
NP_CIF <- function( failtime, status, tt0){
  
  #CIF 
  CIF <- Cuminc( failtime, status )
  
  # CIF for risk 1 and risk 2 at time tt0.
  CIFs_tt0 <- CIF[which.min(abs(tt0 - replace(CIF[,1], CIF[,1]>tt0, Inf))), c(3,4) ]
  
  return( CIFs_tt0 )
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
QP_IRL <- function(data, X, k.num, tt0.in, eps, lambda, eta0){
  
  ##########################
  # Empty vectors/matrices to store values from each repetition.
  M_store <- rep(0,k.num)
  p = length(eta0)
  eta_store <- matrix(0,ncol=p,nrow=k.num)
  #########################
  
  # Clinicians Value function
  CIF_clin <- NP_CIF( data$time, data$status, tt0 = tt0.in )
  
  V_Y_clin <-  as.numeric( CIF_clin[1] )
  V_Z_clin_tol <- as.numeric( CIF_clin[2] ) - lambda 
  
  mu_clin <- c( V_Y_clin, hinge(V_Z_clin_tol) )
  
  # Initialize decision rule
  lin_rule_coeff <- eta0     #c(0.01 , 0 , -0.5 , 0 , 0 , 0) 
  lin_mod1 <- as.vector(lin_rule_coeff %*% t(X))
  rule1 <- trt_rule(lin_mod1)
  
  #rule1
  print(mean(rule1 == data$A))
  
  # Initial estimated rule
  CIFs_eta <- calc.F(eta0, data)
  
  VY_learner <- CIFs_eta$F1t0(tt0.in)
  VZ_learner <- CIFs_eta$F2t0(tt0.in)
  
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
    
    q_coef = 1
    
    # QP set up.
    m = length(mu_clin)
    Dm = diag(m+1) 
    Dm[nrow(Dm)-2, ncol(Dm)-2] <- Dm[nrow(Dm)-1, ncol(Dm)-1] <- 0.000001
    
    dv = rep(0,m+1)
    
    c_m1 <- c(1,1,0)
    c_m2 <- rbind(c(1,0,0),c(0,1,0))
    
    VY_diff <- -(L_mu_mat[1,2] - L_mu_mat[k+1,2])
    VZ_diff <- -(L_mu_mat[1,3] - L_mu_mat[k+1,3])
    V_diff <- c( VY_diff, VZ_diff, q_coef)
    
    Am <- rbind(c_m1,c_m2, V_diff) 
    
    bv <- c(1, 0.00001, 0.00001, 0)
    
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
    plot(L_mu_mat[,-1],col="blue", pch=19, main="Value Function Estimation",cex=1.2)
    points(x= L_mu_mat[1,2],y= L_mu_mat[1,3],bg="red",pch=22,cex=1.5) 
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
    
    fval<-1
    for(j in 1:nrow(eta_start)){
      fit<-try(optim(par=eta_start[j,],fn=opt.eta,mydata=data,tt0=tt0.in,alp=lambda,M=Mk),silent=TRUE)
      if(!is.character(fit)){
        if(fit$value<fval){
          eta<-fit$par/sqrt(sum(fit$par^2))
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
    
    
    eta_k <- eta
    print(eta_k)
    # Store eta values
    eta_store[k,] <- eta_k
    
    # Get decision rule.
    lin_modk <- as.vector(eta_k %*% t(X))
    #print(lin_modk)
    rulek <- trt_rule(lin_modk)
    
    print( c("rule_k",mean(rulek == data$A)) )
    
    CIFs_eta <- calc.F(eta_k, data)
    
    VY_learner <- CIFs_eta$F1t0(tt0.in)
    VZ_learner <- CIFs_eta$F2t0(tt0.in)
    
    # Update Value functions.
    mu_learner <- c(VY_learner, hinge(VZ_learner - lambda) )
    L_mu_mat <- rbind(L_mu_mat, c(-1,mu_learner) )
    
    #print(sum(rulek == data$A) / dim(data)[1])
    # update iteration step.
    
    if(k >= 5 & abs(qk) <= eps ){break}
    
    k = k+1
    
  }
  
  out_list <- list(M_store=M_store, 
                   eta_store = eta_store, 
                   eta_opt = eta_k,
                   M_est = Mk,
                   eps_k = qk,
                   Value_mat = L_mu_mat,
                   trt_match = sum(rulek == data$A) / dim(data)[1],
                   decision_rule = rulek
  )
  
  return(out_list)
  
}

##################################################################################

# real data application
npar<-8
eta0<-matrix(0,2*npar,npar)
for(i in 1:npar){
  eta0[2*i-1,i]<-1
  eta0[2*i,i]<--1
}


# survival time.
tt0.in <- 365*4

# Looking at cumaltive incidence functions to make sure they work.
CIF_clin <- NP_CIF( mydata$time, mydata$status, tt0 = tt0.in )

CIFs_eta <- calc.F(eta0[1,], mydata)

######################################################################################
# Quadratic Programming for AL-IRL.

# survival time.
tt0.in <- 365

xvars <- cbind( rep(1, dim(mydata)[1]), mydata[ , c( "Age", "Male", "Black", "GC", "Single", "Married", "Urban" ) ] )

# Different etas
AL_HIV <- QP_IRL(data=mydata , X = xvars , k.num=10, tt0.in = tt0.in, eps = 0.00001, lambda=0, eta0 = eta0[2,] )

AL_HIV2 <- QP_IRL(data=mydata , X = xvars , k.num=10, tt0.in = tt0.in, eps = 0.00001, lambda=0, eta0 = eta0[6,] )


#########################################################################################


fval<-1
for(j in 1:nrow(eta0)){
  fit<-try(optim(par=eta0[j,],fn=opt.eta,mydata=mydata,tt0=tt0.in,alp=0,M=1000),silent=TRUE)
  if(!is.character(fit)){
    if(fit$value<fval){
      eta<-fit$par/sqrt(sum(fit$par^2))
      #etarec[l,]<-eta
      fval<-fit$value
      cat(paste("j=",j,"\n"))
      print(c(eta,fval))
    }
  }
}


# Compare with large M.
cat("M:", 1000, "eta:", eta)

M1000_rule <- trt_rule(eta %*% t(xvars))

cat("M:", AL_HIV$M_est , "eta:", AL_HIV$eta_opt)

M1_rule <- trt_rule(AL_HIV$eta_opt %*% t(xvars))

sum(M1000_rule == mydata$A)
sum(M1_rule == mydata$A)


