##################################################################################                        
#  Paper 3: Simulation Study                                                     #
#  Rushani Wijesuriya 31 December 2020                                          #
#                                                                               #
#     SL-JM-wide-DI-with mode (waves 2,3 and 4)                                  #
###################################################################################

rm(list=ls())

#load the libraries
library(lme4) #for fitting lmer
library(dplyr)  ##for data wrangling
library(mitml) #for pooling and analysing the results
library(jomo)  #for JM-1L and JM-2L
library(fastDummies) #for creating dummy indicators for the school clusters
library(modeest)  # for finding the mode
library(xlsx)

args = commandArgs(trailingOnly=TRUE)  #these arguments are passed on from a command line that sends script to a HPC server
# Alternatively, can comment out and set parameters as below

datnum<-as.numeric(args[1])   

#load the simulated data
T1 = list.files(pattern=glob2rx("*.csv"))


#select a subset of the data
s2=seq(10,1000,by=10)
s1=s2-9


temp=T1[s1[datnum]:s2[datnum]]

#set seed
set.seed(25218190)
seed=sample(1:1e8,100,replace=F)[datnum] 

set.seed(seed)


results.est=matrix(NA,nrow=6,ncol=length(temp))
results.sd=matrix(NA,nrow=6,ncol=length(temp))
results.RE=matrix(NA,nrow=3,ncol=length(temp))
results.ICC=matrix(NA,nrow=2,ncol=length(temp))
results.CI=c()

data=lapply(temp,read.csv)

for (i in 1: length(temp)){
  print(i)
  dataL=data[[i]]
  
  #reshape to wide format 
  
  dataw <- reshape(dataL,v.names=c("school","t_rating","prev_sdq","prev_dep"),timevar = "time",idvar="c_id",direction= "wide")
 
  
  ##Set number of imputations 
  M=20
  nburn=1000
  NB=10
  
  ##create a dataframe with variables to be imputed
  myvars <- names(dataw) %in% c("prev_dep.2","prev_dep.3","prev_dep.4","prev_dep.5","prev_dep.6","prev_dep.7","prev_dep.8") 
  dataimp=dataw[myvars]
  
  
  ##generate school cluster indicator (mode of school clus at 2,3,and 4)
  School_clusters=dataw %>% select(school.2,school.3,school.4,school.5,school.6,school.7,school.8)
  
  #find the mode and attach it to CATS_data
  dataw$school.mode=apply(School_clusters,1,mfv1)
  
  
  ##create a dataframe with complete variables
  datacomp= cbind(Intercept=rep(1,nrow(dataw)),
                  dataw[,!names(dataw)%in%c("c_id","prev_dep.2","prev_dep.3","prev_dep.4","prev_dep.5","prev_dep.6",
                                            "prev_dep.7","prev_dep.8","school.2",
                                            "school.3","school.4","school.5","school.6","school.7","school.8",
                                            "school.W1" )])
  
  ##create school dummy indicators
  datacomp_dummy=fastDummies::dummy_cols(datacomp, select_columns =c("school.mode"),
                                         remove_first_dummy = TRUE)
  
  #remove school cluster indicators 
  datacomp_dummy <- datacomp_dummy %>% select(-c("school.mode"))
  
  ##perform imputations without the random effects (SL imputation)
  imp1<-jomo1con(Y=dataimp, X=datacomp_dummy, nimp=M,nburn=nburn,nbetween=NB)
  
  # impCheck<-jomo1con.MCMCchain(Y=dataimp, X=datacomp, nburn=nburn)
  # plot(c(1:nburn),impCheck$collectbeta[1,1,1:nburn],type="l")
  
  
  ##Analysis
  
  mylist=list()
  
  for(m in 1:M)
  {
    datw<-imp1[imp1$Imputation==m,]
    
    
    #1. attach school cluster variables and remove DIs
    datw=cbind(datw,dataw$school.2,dataw$school.3,dataw$school.4,dataw$school.5,dataw$school.6,
               dataw$school.7,dataw$school.8)
    
    datw=datw %>% select(!contains(c("school.mode","prev_sdq")),-Imputation,-Intercept)
    
    datw=datw %>% rename(school.2=`dataw$school.2`,school.3=`dataw$school.3`,school.4=`dataw$school.4`,
                         school.5=`dataw$school.5`,school.6=`dataw$school.6`,school.7=`dataw$school.7`,
                         school.8=`dataw$school.8`)
    
    #reshape to long format 
    
    datL=reshape(datw,varying =list(c("prev_dep.2","prev_dep.3","prev_dep.4","prev_dep.5","prev_dep.6","prev_dep.7","prev_dep.8"),
                                    c("t_rating.2","t_rating.3","t_rating.4","t_rating.5","t_rating.6","t_rating.7","t_rating.8"),
                                    c("school.2","school.3","school.4","school.5","school.6","school.7","school.8")),idvar="c_id", 
                 v.names=c("prev_dep","t_rating","school"), times=c(2,3,4,5,6,7,8),direction= "long")
    
    
    #5. save the dataset in a list
    mylist[[m]]= datL
  }
  
  #fit the analysis of interest on the imputed datasets 
  mods <- lapply(mylist,function(d) {lmer(t_rating~ prev_dep+time+c_age+t_score.W1+c_sex+c_ses+
                                            (1|id)+(1|school), data=d) } )
  
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
write.xlsx(results.est,paste0("JM1LDIwide.mode2_Results_est",datnum,".xlsx"))
write.xlsx(results.sd,paste0("JM1LDIwide.mode2_Results_sd",datnum,".xlsx"))
write.xlsx(results.RE,paste0("JM1LDIwide.mode2_Results_RE",datnum,".xlsx"))
write.xlsx(results.CI,paste0("JM1LDIwide.mode2_Results_CI",datnum,".xlsx"))
write.xlsx(results.ICC,paste0("JM1LDIwide.mode2_Results_ICC",datnum,".xlsx"))






