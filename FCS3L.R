##################################################################################                        
#  Paper 3: Simulation Study                                                     #
#  Rushani Wijesuriya 27 December 2020                                          #
#                                                                               #
#     FCS-3L-ml.lmer                                                            #
###################################################################################

rm(list=ls())

#load the libraries
library(lme4) #for fitting lmer
library(dplyr)  ##for data wrangling
library(mitml) #for pooling and analysing the results
library(mice)
library(miceadds)
library(xlsx)


set.seed(228473)

##load the simulated data 
temp = list.files(pattern=glob2rx("*.csv"))

data=lapply(temp,read.csv)

results.est=matrix(NA,nrow=6,ncol=length(temp))
results.sd=matrix(NA,nrow=6,ncol=length(temp))
results.RE=matrix(NA,nrow=3,ncol=length(temp))
results.ICC=matrix(NA,nrow=2,ncol=length(temp))
results.CI=c()

for (i in 1:length(temp)){
  
  print(i)
  dataL=data[[i]]
  
  #remove unwanted variables
  dataL <- dataL %>% select(!"school.W1")
  
  #######specify the imputation model#######################
  
  # create default predictor matrix and imputation methods
  predMatrix <- make.predictorMatrix(data = dataL)
  impMethod <- make.method(data = dataL)
  
  
  # method for lower-level variables :prev_dep at level 1
  impMethod["prev_dep"] <- "ml.lmer"
  # ... for variables at top level (3)
  #none

  # remove cluster indicator variables from set of predictors
  predMatrix[, c("c_id", "school")] <- 0
  
  
  # For variables imputed with ml.lmer, the setup  requires a few extra arguments.
  # Specifically, we need to specify 
  #(a) the level at which the higher-level variables (level 3 and 2) are measured and
  #(b) the cluster variables that define the clustered structure in the imputation model of y.
  
  
  # specify levels of higher-level variables 
  level <- character(ncol(dataL))## change to id
  names(level) <- colnames(dataL)
  level["c_sex"] <- "c_id"    #for sex
  level["c_age"] <- "c_id"      #for age
  level["t_score.W1"] <- "c_id"    #for teacher rating at wave 1
  level["c_ses"] <- "c_id"         # for SES
  
  # specify cluster indicators (as list)
  cluster <- list()
  cluster[["prev_dep"]] <-  c("c_id","school")  #the exposure
  
  imp <- mice(dataL, method = impMethod, predictorMatrix = predMatrix, maxit = 20, 
              m = 20, levels_id = cluster, variables_levels = level,
              aggregate_automatically = FALSE)
  
  
  mylist=list()
  
  for(m in 1:20)
  {
    datw<-mice::complete(imp,m)
    
    #2. remove unwanted variables
    datw=datw %>% select(!contains("sdq"))
    
    
    #5. save the dataset in a list
    mylist[[m]]= datw
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
write.xlsx(results.est,"FCS3L_Results_est.xlsx")
write.xlsx(results.sd,"FCS3L_Results_sd.xlsx")
write.xlsx(results.RE,"FCS3L_Results_RE.xlsx")
write.xlsx(results.CI,"FCS3L_Results_CI.xlsx")
write.xlsx(results.ICC,"FCS3L_Results_ICC.xlsx")

