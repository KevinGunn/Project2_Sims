setwd( ~/EMR-research/Project 2 data )
sepdat<-read.csv( rl_data_final_cont.csv ,head=TRUE)
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

## Pull out first bloc for each subject.

sepdat1 <- subset(sepdat,bloc==1)

# Outcomes: died_in_hosp, reward

# Treatment: IV Fluid and Vasopressor. Recode into new variable.

iv_vp_trt <- rep(0,dim(sepdat1)[1] )

for( i in 1:dim(sepdat1)[1] ){
  
  if( (sepdat1$iv_input[i] == 0) & (sepdat1$vaso_input[i] == 0) ){ 
      # No trt
      iv_vp_trt[i] = 0
      
    } else if ( (sepdat1$iv_input[i] == 1) & (sepdat1$vaso_input[i] == 0) ){
      # iv = 1 and vaso = 0
      iv_vp_trt[i] = 1
      
    } else if ( (sepdat1$iv_input[i] == 2) & (sepdat1$vaso_input[i] == 0) ){
      # iv = 2 and vaso = 0
      iv_vp_trt[i] = 2
      
    } else if ( (sepdat1$iv_input[i] == 3) & (sepdat1$vaso_input[i] == 0) ){
      # iv = 3 and vaso = 0
      iv_vp_trt[i] = 3
      
    } else if ( (sepdat1$iv_input[i] == 4) & (sepdat1$vaso_input[i] == 0) ){
      # iv = 4 and vaso = 0
      iv_vp_trt[i] = 4
      
    } else if ( (sepdat1$iv_input[i] == 0) & (sepdat1$vaso_input[i] == 1) ){
      # iv = 0 and vaso = 1
      iv_vp_trt[i] = 5 
      
    } else if ( (sepdat1$iv_input[i] == 1) & (sepdat1$vaso_input[i] == 1) ){
      # iv = 1 and vaso = 1
      iv_vp_trt[i] = 6 
      
    } else if ( (sepdat1$iv_input[i] == 2) & (sepdat1$vaso_input[i] == 1) ){
      # iv = 2 and vaso = 1
      iv_vp_trt[i] = 7 
      
    } else if ( (sepdat1$iv_input[i] == 3) & (sepdat1$vaso_input[i] == 1) ){
      # iv = 3 and vaso = 1
      iv_vp_trt[i] = 8  
      
    } else if ( (sepdat1$iv_input[i] == 4) & (sepdat1$vaso_input[i] == 1) ){
      # iv = 4 and vaso = 1
      iv_vp_trt[i] = 9  
      
    } else if ( (sepdat1$iv_input[i] == 0) & (sepdat1$vaso_input[i] == 2) ){
      # iv = 0 and vaso = 2
      iv_vp_trt[i] = 10 
      
    } else if ( (sepdat1$iv_input[i] == 1) & (sepdat1$vaso_input[i] == 2) ){
      # iv = 1 and vaso = 2
      iv_vp_trt[i] = 11   
      
    } else if ( (sepdat1$iv_input[i] == 2) & (sepdat1$vaso_input[i] == 2) ){
      # iv = 2 and vaso = 2
      iv_vp_trt[i] = 12 
      
    } else if ( (sepdat1$iv_input[i] == 3) & (sepdat1$vaso_input[i] == 2) ){
      # iv = 3 and vaso = 2
      iv_vp_trt[i] = 13  
      
    } else if ( (sepdat1$iv_input[i] == 4) & (sepdat1$vaso_input[i] == 2) ){
      # iv = 4 and vaso = 2
      iv_vp_trt[i] = 14  
      
    } else if ( (sepdat1$iv_input[i] == 0) & (sepdat1$vaso_input[i] == 3) ){
      # iv = 0 and vaso = 3
      iv_vp_trt[i] = 15   
      
    } else if ( (sepdat1$iv_input[i] == 1) & (sepdat1$vaso_input[i] == 3) ){
      # iv = 1 and vaso = 3
      iv_vp_trt[i] = 16    
      
    } else if ( (sepdat1$iv_input[i] == 2) & (sepdat1$vaso_input[i] == 3) ){
      # iv = 2 and vaso = 3
      iv_vp_trt[i] = 17
      
    } else if ( (sepdat1$iv_input[i] == 3) & (sepdat1$vaso_input[i] == 3) ){
      # iv = 3 and vaso = 3
      iv_vp_trt[i] = 18
      
    } else if ( (sepdat1$iv_input[i] == 4) & (sepdat1$vaso_input[i] == 3) ){
      # iv = 4 and vaso = 3
      iv_vp_trt[i] = 19
      
    } else if ( (sepdat1$iv_input[i] == 0) & (sepdat1$vaso_input[i] == 4) ){
      # iv = 0 and vaso = 4
      iv_vp_trt[i] = 20 
      
    } else if ( (sepdat1$iv_input[i] == 1) & (sepdat1$vaso_input[i] == 4) ){
      # iv = 1 and vaso = 4
      iv_vp_trt[i] = 21 
      
    } else if ( (sepdat1$iv_input[i] == 2) & (sepdat1$vaso_input[i] == 4) ){
      # iv = 2 and vaso = 4
      iv_vp_trt[i] = 22
      
    } else if ( (sepdat1$iv_input[i] == 3) & (sepdat1$vaso_input[i] == 4) ){
      # iv = 3 and vaso = 4
      iv_vp_trt[i] = 23
      
    } else if ( (sepdat1$iv_input[i] == 4) & (sepdat1$vaso_input[i] == 4) ){
      # iv = 4 and vaso = 4
      iv_vp_trt[i] = 24
    
    } 
}    
 

sepdat_out <- data.frame(sepdat1, iv_vp_trt)

hist(sepdat_out$iv_vp_trt)

table(sepdat_out$iv_vp_trt)

sum(table(sepdat_out$iv_vp_trt)[2:5])
#6661
sum(table(sepdat_out$iv_vp_trt)[c(6,11,16,21)])
#173

dim(sepdat_out)

######################################################################################################

# Extract subjects who only received IV_fluids or no treatment.
sepdat_iv <- sepdat[ which(sepdat$vaso_input==0), ] 
  
  
######################################################################################################


thedat <- within(sepdat_iv, {  bloc <- factor(bloc)
iv_input <- factor(iv_input,levels=0:4)
id <- factor(icustayid)
})

#   First create the basic plot object
pp <- ggplot(thedat,aes(x=bloc,y=iv_input,group=id))

#Average Dosage Trajectory using all 20955 subjects.
pp + stat_summary(aes(group = 1,color= red ),
geom =  line , fun.y = mean, size = 1.5) +
xlab( Bloc ) + ylab( Dosage ) + theme(legend.position =  none ) +
ggtitle( Average Trajectory for Subjects Receiving only IV Fluid )
#dev.off()
    
  
# Random subset of patients
sepdat_5 <- sepdat_iv[ which( sepdat_iv$icustayid %in% sepdat_iv$icustayid[sample(nrow(sepdat_iv),5)] ) ,  ]

thedat <- within(sepdat_5, {  bloc <- factor(bloc)
iv_input <- factor(iv_input,levels=0:4)
id <- factor(icustayid)
})

#   First create the basic plot object
pp2 <- ggplot( thedat,aes(x=bloc,y=iv_input,group=id, color=id, shape=id ))

pp2 + geom_line() + geom_point(shape=18) + xlab( Bloc ) + ylab( Dosage ) +  
  theme(legend.position =  none )

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
colnames(trt_bloc_mat) <- c( bloc ,  Low Dose ,  High Dose )

blocs <- sort(unique(sepdat_iv_dose$bloc))

for( i in blocs ){
  
  sepdat_iv_dose_bloc <- subset(sepdat_iv_dose, bloc==i)
  bt <- table(sepdat_iv_dose_bloc$iv_dose)
  trt_bloc_mat[i,] <- c( i, bt[1], bt[2] )
  
}


plot(trt_bloc_mat[,2], ylim = c(0,18000), xlab= Blocks , ylab= Number of Patients )
lines(trt_bloc_mat[,2])
points(trt_bloc_mat[,3], col=2)
lines(trt_bloc_mat[,3],col=2)
# Add a legend
legend( topright , legend=c( Low Dose ,  High Dose ),
       col=c( black ,  red ), lty=1, cex=0.8)

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





