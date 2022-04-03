
##################################################################################                        
#  Paper 3: Simulation Study                                                     #
#  Rushani Wijesuriya 27 January 2021                                          #
#                                                                               #
#     FCS-3L_F-ml.lmer                                                            #
###################################################################################

rm(list=ls())

#load the libraries
library(lme4) #for fitting lmer
library(dplyr)  ##for data wrangling
library(mitml) #for pooling and analysing the results
library(mice)
library(miceadds)
library(modeest)  # for finding the mode
library(xlsx)


set.seed(776564)

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
  
  #reshape to wide format 
  dataw <- reshape(dataL,v.names=c("school","t_rating","prev_sdq","prev_dep"),timevar = "time",idvar="c_id",direction= "wide")
  
  
  ##generate school cluster indicator (mode of school clus at 1,2 and 3)
  School_clusters=dataw %>% select(school.2,school.3,school.4,school.5,school.6,school.7,school.8)
  
  #find the mode and attach it to CATS_data
  dataw$school.mode=apply(School_clusters,1,mfv1)
  
  
  #remove school cluster indicators at waves 1,3 and 4
  dataw <- select(dataw,-c(school.W1,school.2,school.3,school.4,school.5,school.6,school.7,school.8))
  
  datL=reshape(dataw,varying =list(c("prev_dep.2","prev_dep.3","prev_dep.4","prev_dep.5","prev_dep.6","prev_dep.7","prev_dep.8"),
                                   c("t_rating.2","t_rating.3","t_rating.4","t_rating.5","t_rating.6","t_rating.7","t_rating.8"),
                                   c("prev_sdq.2","prev_sdq.3","prev_sdq.4","prev_sdq.5","prev_sdq.6","prev_sdq.7","prev_sdq.8")),
               idvar="c_id", 
               v.names=c("prev_dep","t_rating","prev_sdq"), times=c(2,3,4,5,6,7,8),direction= "long")
  
  datL=datL %>% rename(school=school.mode)
  
  #######specify the imputation model#######################
  
  # create default predictor matrix and imputation methods
  predMatrix <- make.predictorMatrix(data = datL)
  impMethod <- make.method(data = datL)
  
  
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
  level <- character(ncol(datL))## change to id
  names(level) <- colnames(datL)
  level["c_sex"] <- "c_id"    #for sex
  level["c_age"] <- "c_id"      #for age
  level["t_score.W1"] <- "c_id"    #for teacher rating at wave 1
  level["c_ses"] <- "c_id"         # for SES
  
  # specify cluster indicators (as list)
  cluster <- list()
  cluster[["prev_dep"]] <-  c("c_id","school")  #the exposure
  
  imp <- mice(datL, method = impMethod, predictorMatrix = predMatrix, maxit = 20, 
              m = 20, levels_id = cluster, variables_levels = level,
              aggregate_automatically = FALSE)
  
  
  mylist=list()
  
  for(m in 1:20)
  {
    datw<-mice::complete(imp,m)
    
    #2. remove unwanted variables
    datw=datw %>% select(!contains("sdq"))
    
    
    #3. reattach the original school variable
    datw$school <- dataL$school
    
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
write.xlsx(results.est,"FCS3Lc_Results_est.xlsx")
write.xlsx(results.sd,"FCS3Lc_Results_sd.xlsx")
write.xlsx(results.RE,"FCS3Lc_Results_RE.xlsx")
write.xlsx(results.CI,"FCS3Lc_Results_CI.xlsx")
write.xlsx(results.ICC,"FCS3Lc_Results_ICC.xlsx")

