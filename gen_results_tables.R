###################################################################################################################
#          Simulation Study-Summarizing results (Regression coefficient and the variance components)
#          summarizes the performance of each approach in terms of bias,RB modSE, Empse and coverage 
#
###################################################################################################################


rm(list = ls())

library(rsimsum)
library(readxl)
library(WriteXLS)


# ##folder names for changing the directory
MD.mech=c("MAR1","MAR2")

ICC=c("High-high","High-low","Low-high","Low-low")

Clus.size=c("Clus40","Clus10")


#the VC paramaters
VC1=c(0.5,0.8,0.5,0.8)
VC2=c(0.35,0.05,0.45,0.15)
VC3=c(0.15,0.15,0.05,0.05)


for (i in 1:length(MD.mech)) {
  
  for(k in 1:length(ICC)){
    
    for(l in 1:length(Clus.size)){
      
    
      #set path to the working directory with results
      file_path <- file.path(MD.mech[[i]],ICC[[k]],Clus.size[[l]],"/Results")

    ##import all the files in the folder
    temp <- list.files(path=file_path,pattern = glob2rx("*.xlsx")) 
    
    ##extract the results for the estimates only
    temp1=temp[grep("est",temp)]
    
    #  od <- matrix(NA,length(temp1),2)
    #  
    # for (i in 1:length(temp1)){od[i,1] <- temp1[[i]]}
    
    results_prevdep=matrix(NA,nrow=length(temp1),ncol=7)
    
    dup=c("ACA_Results_est.xlsx","FCS1LDIwide.234_Results_est.xlsx",
          "JM1LDIwide.second_Results_est.xlsx",
          "FCS1LDIwide.second_Results_est.xlsx","FCS3Lf_Results_est.xlsx",
          "JM1LDIwide.mode2_Results_est.xlsx",
          "FCS1LDIwide.mode2_Results_est.xlsx","FCS3Lc_Results_est.xlsx","FCS2Lwide.234_Results_est.xlsx",
          "JM2LDI.234_Results_est.xlsx",
          "FCS2LDI.234_Results_est.xlsx","FCS3L_Results_est.xlsx")
    
    temp1=temp1[order(match(temp1, dup))]
    

    dep=c()
    
    for (m in 1:length(temp1)){
      data=read_excel(paste0(file_path,"/",temp1[m]))
      print(dim(data))
      data=data[-1]
      
      dep <- c(dep, as.vector(as.numeric(as.matrix(data[1, ]))))
      
      
    }
    
    ###############################STANDARD ERRORS##########################################################
    
    temp2=temp[grep("sd",temp)]
    
    
    dup=c("ACA_Results_sd.xlsx","FCS1LDIwide.234_Results_sd.xlsx",
          "JM1LDIwide.second_Results_sd.xlsx","FCS1LDIwide.second_Results_sd.xlsx",
          "FCS3Lf_Results_sd.xlsx",
          "JM1LDIwide.mode2_Results_sd.xlsx",
          "FCS1LDIwide.mode2_Results_sd.xlsx","FCS3Lc_Results_sd.xlsx","FCS2Lwide.234_Results_sd.xlsx",
          "JM2LDI.234_Results_sd.xlsx",
          "FCS2LDI.234_Results_sd.xlsx","FCS3L_Results_sd.xlsx")
    
    temp2=temp2[order(match(temp2, dup))]
    
    sd.dep=c()
    
    for (m in 1:length(temp2)){
      data=read_excel(paste0(file_path,"/",temp2[m]))
      data=data[-1]
      print(dim(data))
      sd.dep <- c(sd.dep, as.vector(as.numeric(as.matrix(data[1, ]))))
    }
    
    
    Method.list <- c("ACA","FCS-1L-DI-wide",
                     "JM-1L-DI-wide_f","FCS-1L-DI-wide_f","FCS-3L_f",
                     "JM-1L-DI-wide_c","FCS-1L-DI-wide_c","FCS-3L_c",
                     "FCS-2L-wide","JM-2L-DI",
                     "FCS-2L-DI","FCS-3L")
    
    MSIM_dep <- data.frame(rep(1:1000, length(temp1)),rep(Method.list, each=1000),dep,sd.dep)
    colnames(MSIM_dep)=c("dataset","Method","b","se")
    MSIM_dep <- MSIM_dep[order(MSIM_dep$dataset),]
    
    
    S1=simsum(data=MSIM_dep,estvarname = "b", true =(-0.02), se = "se",methodvar ="Method", ref = "ACA")
    
    # od <- c(1,13,2,8,5,12,10,9,6,3,7,4,11)
    od <- c(1,2,4,7,11,9,12,8,5,6,3,10)
    
    
    ##average estimate
    est=cbind(get_data(S1,stats=c("thetamean")),order=od)
    
    ##bias
    bias=cbind(get_data(S1,stats=c("bias")),order=od)
    
    ##relative bias
    RB= cbind((get_data(S1,stats=c("bias"))[2]/(-0.02))*100,
              order=od)
    
    ##emp standard error
    emp=cbind(get_data(S1,stats=c("empse")),order=od)
    
    #std_bias=cbind((get_data(S1,stats=c("bias"))[2]/get_data(S1,stats=c("empse"))[2])*100,
    #order=c(1,9,7,5,3,8,6,4,2))
    ##model based SE  
    mod=cbind(get_data(S1,stats=c("modelse")),order=od)
    
    ##Coverage 
    cov=cbind(get_data(S1,stats=c("cover")),order=od)
    
    ##Mean Suqared error
    #MSE=cbind(get_data(S1,stats=c("mse")),order=c(1,9,7,5,3,8,6,4,2))
    
    prevdeplist=list(est,bias,RB,emp,mod,cov)
    
    sorted=lapply(prevdeplist,function(d){d[order(d$order),]})
    
    results_prevdep[,1]=c("ACA","FCS-1L-DI-wide",
                          "JM-1L-DI-wide_f","FCS-1L-DI-wide_f","FCS-3L-Blimp_f",
                         "JM-1L-DI-wide_c","FCS-1L-DI-wide_c","FCS-3L-Blimp_c",
                          "FCS-2L-wide","JM-2L-DI",
                          "FCS-2L-DI","FCS-3L-ml.lmer")
    
    ##for prev_dep
    results_prevdep[,2]=round(sorted[[1]][,2],6)
    results_prevdep[,3]=round(sorted[[2]][,2],6)
    results_prevdep[,4]=round(sorted[[3]][,1],6)
    results_prevdep[,5]=round(sorted[[4]][,2],6)
    results_prevdep[,6]=round(sorted[[5]][,2],6)
    results_prevdep[,7]=round(sorted[[6]][,2],6)
    #results[6,3:11]=round(sqrt(sorted[[6]][,2]),6)
    
    results_prevdep=as.data.frame(results_prevdep)
    colnames(results_prevdep)=c("Method","Average estimate","Bias","Relative Bias(%)","Emp.SE","Model SE","coverage")
    
    write.csv(results_prevdep,paste0(file_path,"/","main_effect.Perf.csv"))
    
    ####for the variance components 
    
    ##extract the results for the VC estimates only
    temp3=temp[grep("RE",temp)]
    
    results=matrix(NA,nrow=length(temp3),ncol=10)
    
    dup=c("ACA_Results_RE.xlsx","FCS1LDIwide.234_Results_RE.xlsx",
          "JM1LDIwide.second_Results_RE.xlsx","FCS1LDIwide.second_Results_RE.xlsx","FCS3Lf_Results_RE.xlsx",
          "JM1LDIwide.mode2_Results_RE.xlsx",
          "FCS1LDIwide.mode2_Results_RE.xlsx","FCS3Lc_Results_RE.xlsx",
          "FCS2Lwide.234_Results_RE.xlsx",
          "JM2LDI.234_Results_RE.xlsx",
          "FCS2LDI.234_Results_RE.xlsx","FCS3L_Results_RE.xlsx")
    
    temp3=temp3[order(match(temp3, dup))]
    
    L1=c()
    L2=c()
    L3=c()
    
    for (m in 1:length(temp3)){
      data=read_excel(paste0(file_path,"/",temp3[m]),col_types ="numeric")
      data=data[-1]
      L1 <- c(L1, as.vector(as.matrix(as.numeric(data[3,]))))
      L2 <- c(L2, as.vector(as.matrix(as.numeric(data[2,]))))
      L3 <- c(L3, as.vector(as.matrix(as.numeric(data[1,]))))
    }
    
    MSIM_L1 <- data.frame(rep(1:1000, length(temp3)),rep(Method.list,each = 1000),(L1 * L1),
                          rep(1:1000, length(temp3)))
    
    MSIM_L2 <- data.frame(rep(1:1000, length(temp3)),rep(Method.list,each = 1000),(L2 * L2),
                          rep(1:1000, length(temp3)))
    
    MSIM_L3 <- data.frame(rep(1:1000, length(temp3)),rep(Method.list,each = 1000),(L3 * L3),
                          rep(1:1000, length(temp3)))
    
    
    colnames(MSIM_L1) <- colnames(MSIM_L2) <- colnames(MSIM_L3) <- c("dataset", 
                                                                     "Method", "b", "se")
    
    MSIM_L1 <- MSIM_L1[order(MSIM_L1$dataset),]
    MSIM_L2 <- MSIM_L2[order(MSIM_L2$dataset),]
    MSIM_L3 <- MSIM_L3[order(MSIM_L3$dataset),]
    
    
    ##for level 1
    S1=simsum(data=MSIM_L1,estvarname = "b", true =VC1[k], se = "se",methodvar ="Method", ref = "ACA")
    
    bias=cbind(get_data(S1,stats=c("bias")),order=od)
    
    RB=cbind((get_data(S1,stats=c("bias"))[2]/VC1[k])*100,
             order=od)
    
    emp=cbind(get_data(S1,stats=c("empse")),order=od)
    
    l1list=list(bias,RB,emp)
    
    sorted=lapply(l1list,function(d){d[order(d$order),]})
    
    ##for level 2
    S2=simsum(data=MSIM_L2,estvarname = "b", true =VC2[k], se = "se",methodvar ="Method", ref = "ACA")
    
    bias2=cbind(get_data(S2,stats=c("bias")),order=od)
    
    RB2=cbind((get_data(S2,stats=c("bias"))[2]/VC2[k])*100,
              order=od)
    
    emp2=cbind(get_data(S2,stats=c("empse")),order=od)
    

    l2list=list(bias2,RB2,emp2)
    
    sorted2=lapply(l2list,function(d){d[order(d$order),]})
    
    ##for level 3
    S3=simsum(data=MSIM_L3,estvarname = "b", true =VC3[k], se = "se",methodvar ="Method", ref = "ACA")
    
    bias3=cbind(get_data(S3,stats=c("bias")),order=od)
    
    
    RB3=cbind((get_data(S3,stats=c("bias"))[2]/VC3[k])*100,
              order=od)
    
    emp3=cbind(get_data(S3,stats=c("empse")),order=od)
    

    l3list=list(bias3,RB3,emp3)
    
    sorted3=lapply(l3list,function(d){d[order(d$order),]})
    
    ##Mean Suqared error
    #MSE=cbind(get_data(S1,stats=c("mse")),order=c(1,9,7,5,3,8,6,4,2))
    
    
    results[,1]=c("ACA","FCS-1L-DI-wide",
                  "JM-1L-DI-wide_f","FCS-1L-DI-wide_f","FCS-3L-Blimp_f",
                  "JM-1L-DI-wide_c","FCS-1L-DI-wide_c","FCS-3L-Blimp_c",
                  "FCS-2L-wide","JM-2L-DI",
                  "FCS-2L-DI","FCS-3L-ml.lmer")
    
    ##level 3 fill
    results[,2]=round(sorted3[[1]][,2],6)
    results[,3]=round(sorted3[[2]][,1],6)
    results[,4]=round(sorted3[[3]][,2],6)

    
    ##level 2 fill
    results[,5]=round(sorted2[[1]][,2],6)
    results[,6]=round(sorted2[[2]][,1],6)
    results[,7]=round(sorted2[[3]][,2],6)

    ##level 1 fill
    results[,8]=round(sorted[[1]][,2],6)
    results[,9]=round(sorted[[2]][,1],6)
    results[,10]=round(sorted[[3]][,2],6)

    
    results=as.data.frame(results)
    colnames(results)=c("Method","Bias","Realtive Bias(%)","Emp.SE","Bias","Realtive Bias(%)","Emp.SE",
                        "Bias","Realtive Bias(%)","Emp.SE")
    
    
    write.csv(results,paste0(file_path,"/","variance components.csv"))
    }
    
  }
  
}


################################################################################################################
## New simulation scenarios
################################################################################################################




# ##folder names for changing the directory
sim_scenario=c("stand_alone","small_sample","higher_waves")

theta <- c(-0.02,-0.05,-0.02)

ICC_folder=c("High-high","High-low","Low-high","Low-low")

#the VC paramaters
VC1=c(0.5,0.8,0.5,0.8)
VC2=c(0.35,0.05,0.45,0.15)
VC3=c(0.15,0.15,0.05,0.05)

# 
# 
for (i in 1:length(sim_scenario)) {
  
  for(k in 1:length(ICC_folder)){
    
    #set path to the working directory with results
    file_path <- file.path(sim_scenario[[i]],"MAR2", ICC_folder[[k]],"Results")
    temp <- list.files(path=file_path,pattern = glob2rx("*.xlsx")) 
    
    ##extract the results for the estimates only
    temp1=temp[grep("est",temp)]
    
    results_prevdep=matrix(NA,nrow=length(temp1),ncol=7)
    
    dup=c("ACA_Results_est.xlsx","FCS1LDIwide.234_Results_est.xlsx",
          "JM1LDIwide.second_Results_est.xlsx",
          "FCS1LDIwide.second_Results_est.xlsx","FCS3Lf_Results_est.xlsx",
          "JM1LDIwide.mode2_Results_est.xlsx",
          "FCS1LDIwide.mode2_Results_est.xlsx","FCS3Lc_Results_est.xlsx","FCS2Lwide.234_Results_est.xlsx",
          "JM2LDI.234_Results_est.xlsx",
          "FCS2LDI.234_Results_est.xlsx","FCS3L_Results_est.xlsx")
    
    temp1=temp1[order(match(temp1, dup))]
    
    dep=c()
    
    for (m in 1:length(temp1)){
      
      data=read_excel(paste0(file_path,"/",temp1[m]))
      print(dim(data))
      data=data[-1]
      dep <- c(dep, as.vector(as.numeric(as.matrix(data[1, ]))))
    }
    
    
    ###############################STANDARD ERRORS##########################################################
    
    temp2=temp[grep("sd",temp)]
    
    
    dup=c("ACA_Results_sd.xlsx","FCS1LDIwide.234_Results_sd.xlsx",
          "JM1LDIwide.second_Results_sd.xlsx","FCS1LDIwide.second_Results_sd.xlsx",
          "FCS3Lf_Results_sd.xlsx",
          "JM1LDIwide.mode2_Results_sd.xlsx",
          "FCS1LDIwide.mode2_Results_sd.xlsx","FCS3Lc_Results_sd.xlsx","FCS2Lwide.234_Results_sd.xlsx",
          "JM2LDI.234_Results_sd.xlsx",
          "FCS2LDI.234_Results_sd.xlsx","FCS3L_Results_sd.xlsx")
    
    temp2=temp2[order(match(temp2, dup))]
    
    sd.dep=c()
    
    for (m in 1:length(temp2)){
      data=read_excel(paste0(file_path,"/",temp2[m]))
      data=data[-1]
      print(dim(data))
      sd.dep <- c(sd.dep, as.vector(as.numeric(as.matrix(data[1, ]))))
    }
    
    
    Method.list <- c("ACA","FCS-1L-DI-wide",
                     "JM-1L-DI-wide_f","FCS-1L-DI-wide_f","FCS-3L_f",
                     "JM-1L-DI-wide_c","FCS-1L-DI-wide_c","FCS-3L_c",
                     "FCS-2L-wide","JM-2L-DI",
                     "FCS-2L-DI","FCS-3L")
    
    MSIM_dep <- data.frame(rep(1:1000, length(temp1)),rep(Method.list, each=1000),dep,sd.dep)
    colnames(MSIM_dep)=c("dataset","Method","b","se")
    MSIM_dep <- MSIM_dep[order(MSIM_dep$dataset),]
    
    
    S1=simsum(data=MSIM_dep,estvarname = "b", true =theta[i], se = "se",methodvar ="Method", ref = "ACA")
    
    
    
    # od <- c(1,13,2,8,5,12,10,9,6,3,7,4,11)
    od <- c(1,2,4,7,11,9,12,8,5,6,3,10)
    
    
    ##average estimate
    est=cbind(get_data(S1,stats=c("thetamean")),order=od)
    
    ##bias
    bias=cbind(get_data(S1,stats=c("bias")),order=od)
    
    ##relative bias
    RB= cbind((get_data(S1,stats=c("bias"))[2]/(theta[i]))*100,
              order=od)
    
    ##emp standard error
    emp=cbind(get_data(S1,stats=c("empse")),order=od)
    
    #std_bias=cbind((get_data(S1,stats=c("bias"))[2]/get_data(S1,stats=c("empse"))[2])*100,
    #order=c(1,9,7,5,3,8,6,4,2))
    ##model based SE  
    mod=cbind(get_data(S1,stats=c("modelse")),order=od)
    
    ##Coverage 
    cov=cbind(get_data(S1,stats=c("cover")),order=od)
    
    ##Mean Suqared error
    #MSE=cbind(get_data(S1,stats=c("mse")),order=c(1,9,7,5,3,8,6,4,2))
    
    prevdeplist=list(est,bias,RB,emp,mod,cov)
    
    sorted=lapply(prevdeplist,function(d){d[order(d$order),]})
    
    results_prevdep[,1]=c("ACA","FCS-1L-DI-wide",
                          "JM-1L-DI-wide_f","FCS-1L-DI-wide_f","FCS-3L_f",
                          "JM-1L-DI-wide_c","FCS-1L-DI-wide_c","FCS-3L_c",
                          "FCS-2L-wide","JM-2L-DI",
                          "FCS-2L-DI","FCS-3L")
    
    ##for prev_dep
    results_prevdep[,2]=round(sorted[[1]][,2],6)
    results_prevdep[,3]=round(sorted[[2]][,2],6)
    results_prevdep[,4]=round(sorted[[3]][,1],6)
    results_prevdep[,5]=round(sorted[[4]][,2],6)
    results_prevdep[,6]=round(sorted[[5]][,2],6)
    results_prevdep[,7]=round(sorted[[6]][,2],6)
    #results[6,3:11]=round(sqrt(sorted[[6]][,2]),6)
    
    results_prevdep=as.data.frame(results_prevdep)
    colnames(results_prevdep)=c("Method","Average estimate","Bias","Relative Bias(%)","Emp.SE","Model SE","coverage")
    
    write.csv(results_prevdep,paste0(file_path,"/","main_effect.Perf.csv"))
    
    ####for the variance components 
    
    ##extract the results for the VC estimates only
    temp3=temp[grep("RE",temp)]
    
    results=matrix(NA,nrow=length(temp3),ncol=10)
    
    dup=c("ACA_Results_RE.xlsx","FCS1LDIwide.234_Results_RE.xlsx",
          "JM1LDIwide.second_Results_RE.xlsx","FCS1LDIwide.second_Results_RE.xlsx","FCS3Lf_Results_RE.xlsx",
          "JM1LDIwide.mode2_Results_RE.xlsx",
          "FCS1LDIwide.mode2_Results_RE.xlsx","FCS3Lc_Results_RE.xlsx",
          "FCS2Lwide.234_Results_RE.xlsx",
          "JM2LDI.234_Results_RE.xlsx",
          "FCS2LDI.234_Results_RE.xlsx","FCS3L_Results_RE.xlsx")
    
    temp3=temp3[order(match(temp3, dup))]
    
    
    L1=c()
    L2=c()
    L3=c()
    
    for (m in 1:length(temp3)){
      data=read_excel(paste0(file_path,"/",temp3[m]),col_types ="numeric")
      data=data[-1]
      L1 <- c(L1, as.vector(as.matrix(as.numeric(data[3,]))))
      L2 <- c(L2, as.vector(as.matrix(as.numeric(data[2,]))))
      L3 <- c(L3, as.vector(as.matrix(as.numeric(data[1,]))))
      
    }
    
    
    MSIM_L1 <- data.frame(rep(1:1000, length(temp3)),rep(Method.list,each = 1000),(L1 * L1),
                          rep(1:1000, length(temp3)))
    
    MSIM_L2 <- data.frame(rep(1:1000, length(temp3)),rep(Method.list,each = 1000),(L2 * L2),
                          rep(1:1000, length(temp3)))
    
    MSIM_L3 <- data.frame(rep(1:1000, length(temp3)),rep(Method.list,each = 1000),(L3 * L3),
                          rep(1:1000, length(temp3)))
    
    colnames(MSIM_L1) <- colnames(MSIM_L2) <- colnames(MSIM_L3) <- c("dataset", 
                                                                     "Method", "b", "se")
    
    
    MSIM_L1 <- MSIM_L1[order(MSIM_L1$dataset),]
    MSIM_L2 <- MSIM_L2[order(MSIM_L2$dataset),]
    MSIM_L3 <- MSIM_L3[order(MSIM_L3$dataset),]
    
    ##for level 1
    S1=simsum(data=MSIM_L1,estvarname = "b", true =VC1[k], se = "se",methodvar ="Method", ref = "ACA")
    
    bias=cbind(get_data(S1,stats=c("bias")),order=od)
    
    RB=cbind((get_data(S1,stats=c("bias"))[2]/VC1[k])*100,
             order=od)
    
    emp=cbind(get_data(S1,stats=c("empse")),order=od)
    
    l1list=list(bias,RB,emp)
    
    sorted=lapply(l1list,function(d){d[order(d$order),]})
    
    ##for level 2
    S2=simsum(data=MSIM_L2,estvarname = "b", true =VC2[k], se = "se",methodvar ="Method", ref = "ACA")
    
    bias2=cbind(get_data(S2,stats=c("bias")),order=od)
    
    RB2=cbind((get_data(S2,stats=c("bias"))[2]/VC2[k])*100,
              order=od)
    
    emp2=cbind(get_data(S2,stats=c("empse")),order=od)
    
    
    l2list=list(bias2,RB2,emp2)
    
    sorted2=lapply(l2list,function(d){d[order(d$order),]})
    
    ##for level 3
    S3=simsum(data=MSIM_L3,estvarname = "b", true =VC3[k], se = "se",methodvar ="Method", ref = "ACA")
    
    bias3=cbind(get_data(S3,stats=c("bias")),order=od)
    
    
    RB3=cbind((get_data(S3,stats=c("bias"))[2]/VC3[k])*100,
              order=od)
    
    emp3=cbind(get_data(S3,stats=c("empse")),order=od)
    
    
    l3list=list(bias3,RB3,emp3)
    
    sorted3=lapply(l3list,function(d){d[order(d$order),]})
    
    ##Mean Suqared error
    #MSE=cbind(get_data(S1,stats=c("mse")),order=c(1,9,7,5,3,8,6,4,2))
    
    
    results[,1]=c("ACA","FCS-1L-DI-wide",
                  "JM-1L-DI-wide_f","FCS-1L-DI-wide_f","FCS-3L_f",
                  "JM-1L-DI-wide_c","FCS-1L-DI-wide_c","FCS-3L_c",
                  "FCS-2L-wide","JM-2L-DI",
                  "FCS-2L-DI","FCS-3L")
    
    ##level 3 fill
    results[,2]=round(sorted3[[1]][,2],6)
    results[,3]=round(sorted3[[2]][,1],6)
    results[,4]=round(sorted3[[3]][,2],6)
    
    
    ##level 2 fill
    results[,5]=round(sorted2[[1]][,2],6)
    results[,6]=round(sorted2[[2]][,1],6)
    results[,7]=round(sorted2[[3]][,2],6)
    
    ##level 1 fill
    results[,8]=round(sorted[[1]][,2],6)
    results[,9]=round(sorted[[2]][,1],6)
    results[,10]=round(sorted[[3]][,2],6)
    
    
    results=as.data.frame(results)
    colnames(results)=c("Method","Bias","Realtive Bias(%)","Emp.SE","Bias","Realtive Bias(%)","Emp.SE",
                        "Bias","Realtive Bias(%)","Emp.SE")
    
    
    write.csv(results,paste0(file_path,"/","variance components.csv"))
  }
  
}


