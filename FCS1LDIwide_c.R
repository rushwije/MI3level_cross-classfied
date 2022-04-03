##################################################################################                        
#  Paper 3: Simulation Study                                                     #
#  Rushani Wijesuriya 31 December 2020                                          #
#                                                                               #
#     SL-FCS-wide-DI-with mode of waves 2,3 and 4                                #
###################################################################################

rm(list=ls())

#load the libraries
library(lme4) #for fitting lmer
library(dplyr)  ##for data wrangling
library(mitml) #for pooling and analysing the results
library(jomo)  #for JM-1L and JM-2L
library(fastDummies) #for creating dummy indicators for the school clusters
library(modeest)  # for finding the mode
library(mice)
library(xlsx)

set.seed(825374)

##load the simulated data 
temp = list.files(pattern=glob2rx("*.csv"))


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

  ##generate school cluster indicator (mode of school clus at 1,2 and 3)
  School_clusters=dataw %>%  select(school.2,school.3,school.4,school.5,school.6,school.7,school.8)
  
  #find the mode and attach it to CATS_data
  dataw$school.mode=apply(School_clusters,1,mfv1)
  
  ##create school dummy indicators
  dataw$school.mode=as.factor(dataw$school.mode)
  
  
  #specify the predictor matrix
  pred=mice::make.predictorMatrix(dataw)
  pred[!names(dataw)%in% c("prev_dep.2","prev_dep.3","prev_dep.4","prev_dep.5","prev_dep.6","prev_dep.7","prev_dep.8"),]<-0
  pred[,"c_id"] <- 0
  
  #IMPORTANT: here the school cluster indicator at the current wave is used (dep.1 <- cluster 1,dep.2 <- cluster 2 etc)
  pred[c("prev_dep.2","prev_dep.3","prev_dep.4","prev_dep.5","prev_dep.6","prev_dep.7","prev_dep.8"),
       c("school.W1","school.2","school.3","school.4","school.5","school.6","school.7","school.8")]<-0
  
  #specify the methods 
  meth=make.method(dataw)
  meth[c( "prev_dep.2","prev_dep.3","prev_dep.4","prev_dep.5","prev_dep.6","prev_dep.7","prev_dep.8")]="norm"
  
  
  ##Set number of imputations and number of burn-in iterations
  M<-20
  
  
  ##Perform imputations
  imp2=mice(data=dataw,predictorMatrix = pred,m=M,method=meth,maxit=10)
  
  
  ##Analysis
  
  mylist=list()
  
  for(m in 1:M)
  {
    datw<-mice::complete(imp2,m)
    
    #remove unwanted variables
    datw <- datw %>% select(!contains("prev_sdq")& !school.mode)
    
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
                                            (1|c_id)+(1|school), data=d) } )
  
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
write.xlsx(results.est,"FCS1LDIwide.mode2_Results_est.xlsx")
write.xlsx(results.sd,"FCS1LDIwide.mode2_Results_sd.xlsx")
write.xlsx(results.RE,"FCS1LDIwide.mode2_Results_RE.xlsx")
write.xlsx(results.CI,"FCS1LDIwide.mode2_Results_CI.xlsx")
write.xlsx(results.ICC,"FCS1LDIwide.mode2_Results_ICC.xlsx")




