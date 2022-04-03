##################################################################################                        
#  Paper 3: Simulation Study                                                     #
#  Rushani Wijesuriya 25 December 2020                                          #
#                                                                               #
#     ACA                                                                       #
###################################################################################

rm(list=ls())

#load the libraries
library(lme4) #for fitting lmer
library(dplyr)  ##for data wrangling
library(mitml) #for pooling and analysing the results
library(jomo)  #for JM-1L and JM-2L
library(fastDummies) #for creating dummy indicators for the school clusters
library(xlsx)

set.seed(3761911)

##load the simulated data 
temp = list.files(pattern=glob2rx("*.csv"))

data=lapply(temp,read.csv)

comp_results.est=matrix(NA,nrow=6,ncol=length(temp))
comp_results.sd=matrix(NA,nrow=6,ncol=length(temp))
comp_results.RE=matrix(NA,nrow=3,ncol=length(temp))
comp_results.ICC=matrix(NA,nrow=2,ncol=length(temp))
comp_results.CI=c()



for (i in 1:length(temp)){
  
 
  
  dataL=data[[i]]
  fit<-lmer(t_rating~ prev_dep+time+c_age+t_score.W1+c_sex+c_ses+
                      (1|c_id)+(1|school), data=dataL) 
  x=as.data.frame(VarCorr(fit),comp=c("Variance"))
  y=confint(fit)
  comp_results.est[,i]=coef(summary(fit))[2:7,1]
  comp_results.sd[,i]<-coef(summary(fit))[2:7,2]
  comp_results.RE[1,i]=x[2,5]
  comp_results.RE[2,i]=x[1,5]
  comp_results.RE[3,i]=x[3,5]
  comp_results.ICC[1,i]=x[2,4]/(x[2,4]+x[1,4]+x[3,4])
  comp_results.ICC[2,i]=(x[2,4]+x[1,4])/(x[2,4]+x[1,4]+x[3,4])
  comp_results.CI=cbind(comp_results.CI,y[5:10,1:2])
  
}

rownames(comp_results.est)=c("prev_dep","time","c_age","tratingW1","c_gender","c_ses")
colnames(comp_results.est)=c(1:length(temp))

rownames(comp_results.sd)=c("prev_dep","time","c_age","tratingW1","c_gender","c_ses")
colnames(comp_results.sd)=c(1:length(temp))

rownames(comp_results.RE)=c("level 3","level 2","level 1")
colnames(comp_results.RE)=c(1:length(temp))

rownames(comp_results.ICC)=c("level 3","level 2")
colnames(comp_results.ICC)=c(1:length(temp))

rownames(comp_results.CI)=c("prev_dep","time","c_age","tratingW1","c_gender","c_ses")

colnames(comp_results.CI)=c(rep(c(1:length(temp)), each=2))



write.xlsx(comp_results.est,"ACA_Results_est.xlsx")
write.xlsx(comp_results.sd,"ACA_Results_sd.xlsx")
write.xlsx(comp_results.RE,"ACA_Results_RE.xlsx")
write.xlsx(comp_results.CI,"ACA_Results_CI.xlsx")
write.xlsx(comp_results.ICC,"ACA_Results_ICC.xlsx")



  
  
