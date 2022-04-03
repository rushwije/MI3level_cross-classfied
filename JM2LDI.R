##################################################################################                        
#  Paper 3: Simulation Study                                                     #
#  Rushani Wijesuriya 28 December 2020                                          #
#                                                                               #
#     JM-2L-DI WITH WAVE 2,3 and 4                                               #
###################################################################################


rm(list=ls())

#load the libraries
library(lme4) #for fitting lmer
library(dplyr)  ##for data wrangling
library(mitml) #for pooling and analysing the results
library(mice)
library(xlsx)
library(jomo)

args = commandArgs(trailingOnly=TRUE)  #these arguments are passed on from a command line that sends script to a HPC server
# Alternatively, can comment out and set parameters as below

datnum<-as.numeric(args[1])   

#load the simulated data
T1 = list.files(pattern=glob2rx("*.csv"))


#select a subset of the data
s2=seq(100,1000,by=100)
s1=s2-99

temp=T1[s1[datnum]:s2[datnum]]

set.seed(88734)
seed=sample(1:1e8,10,replace=F)[datnum] 

set.seed(seed)

data=lapply(temp,read.csv)

results.est=matrix(NA,nrow=6,ncol=length(temp))
results.sd=matrix(NA,nrow=6,ncol=length(temp))
results.RE=matrix(NA,nrow=3,ncol=length(temp))
results.ICC=matrix(NA,nrow=2,ncol=length(temp))
results.CI=c()


for (i in 1:length(temp)){
  
  print(i)
  dataL=data[[i]]
  
  ##Set number of imputations 
  M=20
  nburn=1000
  NB=100
  
  #################JM-2L-DI wide waves 2, 3 and 4################
  
  ##create a dataframe with variables to be imputed
  myvars <- names(dataL) %in% c("prev_dep") 
  dataimp=dataL[myvars]
  
  ##create a dataframe with complete variables
  datacomp=cbind(Intercept=rep(1,nrow(dataL)),
                 dataL[,!names(dataL)%in%c("prev_dep","c_id","school.W1")])
  
  ##create school dummy indicators for waves, 2,3 and 4
  datacomp_dummy=fastDummies::dummy_cols(datacomp, select_columns =c("school"),
                                         remove_first_dummy = TRUE)
  
  #remove school cluster indicator 
  datacomp_dummy <- datacomp_dummy %>% select(!school)
  
  
  #create a dataframe with random effects
  datLcompRE<-cbind(Intercept=rep(1,nrow(dataL)))

  ##perform imputations with the random effects (2L imputation)
  imp1<-jomo(Y=dataimp, X=datacomp_dummy, clus=dataL$c_id,Z=datLcompRE,nimp=M,nburn=nburn,nbetween=NB)
  
  # impCheck<-jomo1con.MCMCchain(Y=dataimp, X=datacomp, nburn=nburn)
  # plot(c(1:nburn),impCheck$collectbeta[1,1,1:nburn],type="l")
  
  ##Analysis
  
  mylist=list()
  
  for(m in 1:M)
  {
    datL<-imp1[imp1$Imputation==m,]
    
    datL <- datL[!duplicated(as.list(datL))]
    
    #1. remove unwanted variables
    datL=datL %>% select(!contains(c("school_","prev_sdq")),-Intercept,-id)
    
    #reattach the school varible
    datL$school <- dataL$school
    
    #5. save the dataset in a list
    mylist[[m]]= datL
  }
  
  #fit the analysis of interest on the imputed datasets 
  mods <- lapply(mylist,function(d) {lmer(t_rating~ prev_dep+time+c_age+t_score.W1+c_sex+c_ses+
                                            (1|clus)+(1|school), data=d) } )
  
  MI_est=testEstimates(mods, var.comp=TRUE)
  
  #Compute CIs
  CI=matrix(NA,6,2)
  Conf=confint(MI_est)
  CI[,1]=as.vector(Conf[2:7,1])
  CI[,2]=as.vector(Conf[2:7,2])
  results.CI=cbind(results.CI,CI)
  
  
  #store the estimates
  results.est[,i]=MI_est$estimates[2:7,1]
  results.sd[,i]=MI_est$estimates[2:7,2]
  results.RE[1,i]=sqrt(MI_est$var.comp[2,1])
  results.RE[2,i]=sqrt(MI_est$var.comp[1,1])
  results.RE[3,i]=sqrt(MI_est$var.comp[3,1])
  
  
  results.ICC[1,i]=MI_est$var.comp[2,1]/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  results.ICC[2,i]=(MI_est$var.comp[2,1]+MI_est$var.comp[1,1])/(MI_est$var.comp[2,1]+MI_est$var.comp[1,1]+MI_est$var.comp[3,1])
  
}

rownames(results.est)=c("prev_dep","time","c_age","tratingW1","c_gender","c_ses")
colnames(results.est)=c(seq(1:length(temp)))

rownames(results.sd)=c("prev_dep","time","c_age","tratingW1","c_gender","c_ses")
colnames(results.sd)=c(seq(1:length(temp)))

rownames(results.RE)=c("level 3","level 2","level 1")
colnames(results.RE)=c(seq(1:length(temp)))

rownames(results.ICC)=c("level 3","level 2")
colnames(results.ICC)=c(seq(1:length(temp)))

rownames(results.CI)=c("prev_dep","time","c_age","tratingW1","c_gender","c_ses")
colnames(results.CI)=c(rep(1:length(temp), each=2))

##save the results
write.xlsx(results.est,paste0("JM2LDI.234_Results_est",datnum,".xlsx"))
write.xlsx(results.sd,paste0("JM2LDI.234_Results_sd",datnum,".xlsx"))
write.xlsx(results.RE,paste0("JM2LDI.234_Results_RE",datnum,".xlsx"))
write.xlsx(results.CI,paste0("JM2LDI.234_Results_CI",datnum,".xlsx"))
write.xlsx(results.ICC,paste0("JM2LDI.234_Results_ICC",datnum,".xlsx"))

