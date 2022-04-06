##################################################################################                                                                      #
#  Simulation of complete and missing data for the simulation scenarios           #
#  evaluted in the paper :                                                        #
#  "Multiple imputation approaches for handling incomplete three-level data       #
#   with time-varying cluster-memberships"                                        #
#   Author: Rushani Wijesuriya                                                    #
#   08th of December 2020                                                         #
##################################################################################

rm(list=ls())

library(ReIns) #for generating random numbers from log-normal distribution 
library(splitstackshape) #to exapand the rows
library(boot)#for inverting the logit function 
library(dplyr)
library(DataCombine) #for the slide function 


#Function for generating the data 


#seed <-  random state
#simno <- number of replications
#I <- number of higher-level clsuters at wave 1
#J <- average number of students per school
#mdm <- missing data mechanism
#ICCL2 <- Intra cluster correlation at level 2
#ICCL3 <- Intra cluster correlation at level 3
#nw <- no of new higher-level clusters being added at each wave
#waveno <- number of waves of data collection

simul <- function(seed,simno,I,J,mdm,ICCL3,ICCL2,nw,waveno){

set.seed(seed)


#############parameter values#############

#sample size
N <- I*J

#child's age
a <- 7
b <- 10

#parameters for log school frequencies 
log.mu <- 3.27
log.sigma <- 0.57
log.min <- 2.10
log.max <- 4.20

##child's gender
lambda <- 0.5

##chid's teacher rating at wave 1
eta0 <- 3.00
eta1 <- 0.20
eta2 <- 0.10
eta3 <- 0.02

#depressive symptom scores at waves 1,2 and 3 (exposure)
delta0 <- 6.0
delta1 <- 0.01
delta2 <- -0.04
delta3 <- -0.20
delta4 <- -0.10
delta5 <- -0.30

#teacher ratings waves 2,3 and 4 (outcome)
gamma0 <- 2.60
if (N==1200){gamma1 <- -0.02} else if (N==300){gamma1 <- -0.05}
gamma2 <- -0.20
gamma3 <- 0.02
gamma4 <- 0.70
gamma5 <- 0.01
gamma6 <- 0.01

#SDQ values at waves 1,2 and 3 (auxiliary)
beta0 <- 16.0
beta1 <- 0.40
beta2 <- 0.05

##parameters for the random effects

#exposure-depressive symptoms
sd_depL3 <- 0.25
sd_depL2 <- 0.80
sd_depL1 <- 1.50

#outcome-Teacher ratings

if(ICCL3=="High" & ICCL2=="high"){
  sd_tratingL3 <- sqrt(0.15)
  sd_tratingL2 <- sqrt(0.35)
  sd_tratingL1 <-sqrt(0.5)} else if (ICCL3=="High" & ICCL2=="low") {
  sd_tratingL3 <- sqrt(0.15)
  sd_tratingL2 <- sqrt(0.05)
  sd_tratingL1 <-sqrt(0.8)} else if (ICCL3=="Low" & ICCL2=="high") {
  sd_tratingL3 <- sqrt(0.05)
  sd_tratingL2 <- sqrt(0.04)
  sd_tratingL1 <-sqrt(0.5)} else if (ICCL3=="Low" & ICCL2=="low") {
  sd_tratingL3 <- sqrt(0.05)
  sd_tratingL2 <- sqrt(0.15)
  sd_tratingL1 <-sqrt(0.8)} 

#Auxiliary-SDQ
sd_sdqL3 <- 0.8
sd_sdqL2 <- 4.0
sd_sdqL1 <- 2.8

#porportion of students who move (at wave 2,3 and 4)
prop <- 0.05


if(mdm=="MAR1"){
  iota01 <- -18.5
  iota02 <- -19.0
  iota03 <- -21.0
  iota1 <- 1.4
  iota2 <- 1.1 } else if (mdm=="MAR2" & waveno==4){
    iota01 <- -18.5
    iota02 <- -19.0
    iota03 <- -21.0
    iota1 <- 1.4
    iota2 <- 1.1
  }else if (mdm=="MAR2" & waveno==8){
    iota0 <- -16.5
    iota1 <- 1.4
    iota2 <- 1.1
  }



########### Complete data generation #####################


for(i in 1:simno){
  
  ##generate the school clusters at wave 1
  D <- data.frame(school.1=c(1:I))
  
  #generate school level RE for wave 1
  D$L3_RE_dep.1 <- rnorm(I,0,sd_depL3)
  D$L3_RE_trating.1 <- rnorm(I,0,sd_tratingL3)
  D$L3_RE_sdq.1 <- rnorm(I,0,sd_sdqL3)
  
  
  if (nw!=0){
    
  
  #Generate the class sizes
  D$freq <-round(exp(rtlnorm(I,log.mu,log.sigma,endpoint = log.max)))
  
  #obtain the total number of students
  total <- sum(D$freq)
  
  #recale to 1200
  D$class_size <- round(D$freq*(N/total))
  
  #add/remove students from the last cluster to total the no of students to 1200
  
  if ((sum(D$class_size)-N)<0){D$class_size[I] <- D$class_size[I]+(N-sum(D$class_size))} 
  if((sum(D$class_size)-N)>0){D$class_size[I] <- D$class_size[I]-(sum(D$class_size)-N)} 
  
  } else if (nw==0){
    
    N <- I*J
    D$class_size <- J
    
  }
  #replicate the rows for adding students
  D <- expandRows(D, "class_size")
  
  #generate the indivdiuals
  D$c_id <- c(1:N)
  
  D <- D[order(D$school.1),]
  
  D$freq=NULL
  
  #Generate wave 2 school cluster memberships (add 50 schools and assign 5% of students to these schools)
  
  # function to randomly split data 
  random.split <- function(p){
    picked <- sample(seq_len(nrow(D)),size = floor(p*nrow(D)))
    obj <- list(change <- D[picked,], constant <- D[-picked,])
    return(obj)
  }
  
  if (nw!=0){
  
  for (m in 2:waveno){
    
  W <- random.split(prop)
  change <- W[[1]]
  
  #assign the students who moved to 50 new schools
  school <- rep((I+1):(I+nw), length=nrow(change))
  
  #generate school level RE for wave 2 for those who moved
  L3_RE_dep <- rep(rnorm(nw,0,sd_depL3), length=nrow(change))
  L3_RE_trating <- rep(rnorm(nw,0,sd_tratingL3),length=nrow(change))
  L3_RE_sdq <- rep(rnorm(nw,0,sd_sdqL3),length=nrow(change))
  
  DF <- data.frame(school,L3_RE_dep,L3_RE_trating,L3_RE_sdq)

  names(DF) <- c(paste0('school.',m),paste0('L3_RE_dep.',m),paste0('L3_RE_trating.',m),paste0('L3_RE_sdq.',m))
  
  change <- cbind(change, DF)
  
  
  #generate(repeat wave 1)wave 2 information for those who didnt move
  constant <- W[[2]]
  new <- select(constant, ends_with(paste0(m-1)))
  names(new) <- c(paste0('school.',m),paste0('L3_RE_dep.',m),paste0('L3_RE_trating.',m),paste0('L3_RE_sdq.',m))
  
  constant <- cbind(constant, new)
  
  #merge the data 
  D <- rbind(constant,change)
  }
  
  } else if (nw==0){
    
    for (m in 2:waveno){
      
      W <- random.split(prop)
      change <- W[[1]]
      
      #assign the students moving to different schools
      

      new_school_list <- c()
      for (j in 1:nrow(change)){
        
        new_school <- sample(1:10,1)
        
        pre_school <- change[,names(change)%in%c(paste0("school.",m-1))] 
                                    
        while(pre_school[j]==new_school){
          new_school <- sample(1:10,1)
        }
        
        new_school_list[j] <- new_school
      }
      
      change <- cbind(change,new_school_list)
      
      names(change)[length(names(change))] <- paste0("school.",m)
      #generate school level RE for wave 2 for those who moved
      
      constant_schools <- D[,1:4]
      constant_schools <- constant_schools %>% distinct()
      
      names(constant_schools) <- c(paste0('school.',m),paste0('L3_RE_dep.',m),paste0('L3_RE_trating.',m),paste0('L3_RE_sdq.',m))
  
    
      change <- left_join(change,constant_schools,by=paste0("school.",m))
      
      #generate(repeat wave 1)wave 2 information for those who didnt move
      constant <- W[[2]]
      new <- select(constant, ends_with(paste0(m-1)))
      names(new) <- c(paste0('school.',m),paste0('L3_RE_dep.',m),paste0('L3_RE_trating.',m),paste0('L3_RE_sdq.',m))
      
      constant <- cbind(constant, new)
      
      #merge the data 
      D <- rbind(constant,change)
      
    }
    
  }
 
  #rearrange columns
  
  if (waveno==4){D <- D[,c( "c_id","school.1", "L3_RE_dep.1","L3_RE_trating.1","L3_RE_sdq.1","school.2","L3_RE_dep.2",
                            "L3_RE_trating.2","L3_RE_sdq.2","school.3","L3_RE_dep.3","L3_RE_trating.3","L3_RE_sdq.3","school.4",        
                            "L3_RE_dep.4","L3_RE_trating.4","L3_RE_sdq.4")]
  } else if (waveno==8){D <- D[,c("c_id","school.1", "L3_RE_dep.1","L3_RE_trating.1","L3_RE_sdq.1","school.2","L3_RE_dep.2",
                                  "L3_RE_trating.2","L3_RE_sdq.2","school.3","L3_RE_dep.3","L3_RE_trating.3","L3_RE_sdq.3","school.4",        
                                  "L3_RE_dep.4","L3_RE_trating.4","L3_RE_sdq.4","school.5", "L3_RE_dep.5","L3_RE_trating.5","L3_RE_sdq.5",
                                  "school.6", "L3_RE_dep.6","L3_RE_trating.6","L3_RE_sdq.6","school.7", "L3_RE_dep.7","L3_RE_trating.7","L3_RE_sdq.7",
                                  "school.8", "L3_RE_dep.8","L3_RE_trating.8","L3_RE_sdq.8")]}
  
  D <- D[order(D$c_id),]
  
  #child's age at wave 1
  D$c_age <- runif((I*J),a,b)
  
  #Simulate sex (M,F) groups (males=1, females=0)
  D$uran=runif((I*J),0,1)
  
  D$c_sex=ifelse(D$uran<=lambda,1,0)
  
  ##simulate SES values
  D$c_ses <- rnorm(N,0,1)
  
  #Simulate teacher ratings at wave 1
  e_teacherw1 <- rnorm(N,0,1)
  D$t_score.W1=eta0+eta1*D$c_sex+eta2*D$c_ses+eta3*D$c_age+e_teacherw1
  
  D$uran=NULL
  
  #generate individual level REs
  D$L2_RE_dep <- rnorm(N,0,sd_depL2)
  D$L2_RE_trating <- rnorm(N,0,sd_tratingL2)
  D$L2_RE_sdq <- rnorm(N,0,sd_sdqL2)
  
  #reshape to long format
  
  if (waveno==4){ D_long=reshape(D,varying =list(c("school.1","school.2","school.3","school.4"),
                                                 c("L3_RE_dep.1","L3_RE_dep.2","L3_RE_dep.3","L3_RE_dep.4"),
                                                 c("L3_RE_trating.1","L3_RE_trating.2","L3_RE_trating.3","L3_RE_trating.4"),
                                                 c("L3_RE_sdq.1","L3_RE_sdq.2","L3_RE_sdq.3","L3_RE_sdq.4")),idvar="c_id", 
                                 v.names=c("school","L3_RE_dep","L3_RE_trating","L3_RE_sdq"), times=c(1,2,3,4),direction= "long")
  } else if (waveno==8){ D_long=reshape(D,varying =list(c("school.1","school.2","school.3","school.4","school.5","school.6","school.7","school.8"),
                                                        c("L3_RE_dep.1","L3_RE_dep.2","L3_RE_dep.3","L3_RE_dep.4","L3_RE_dep.5","L3_RE_dep.6","L3_RE_dep.7","L3_RE_dep.8"),
                                                        c("L3_RE_trating.1","L3_RE_trating.2","L3_RE_trating.3","L3_RE_trating.4","L3_RE_trating.5","L3_RE_trating.6","L3_RE_trating.7","L3_RE_trating.8"),
                                                        c("L3_RE_sdq.1","L3_RE_sdq.2","L3_RE_sdq.3","L3_RE_sdq.4","L3_RE_sdq.5","L3_RE_sdq.6","L3_RE_sdq.7","L3_RE_sdq.8")),idvar="c_id", 
                                        v.names=c("school","L3_RE_dep","L3_RE_trating","L3_RE_sdq"), times=c(1,2,3,4,5,6,7,8),direction= "long")
  } 
 
  
  D_long <- D_long[order(D_long$c_id),]
  
  #D_long <- D_long %>%filter(time!=1) 
  
  #generate the exposure (depressive symptoms at waves 1,2 and 3)
  e_dep <- rnorm(dim(D_long)[1],0,sd_depL1)
  D_long$dep <- delta0+delta1*D_long$c_age+delta2*D_long$c_sex+delta3*D_long$t_score.W1+
    delta4*D_long$c_ses+delta5*D_long$time+D_long$L3_RE_dep+D_long$L2_RE_dep+e_dep
  
  #generate previous wave depression
  D_long<- slide(D_long, Var = "dep", GroupVar = "c_id",
                 slideBy = -1)
  
  colnames(D_long)[colnames(D_long)=="dep-1"] <- "prev_dep"
  
  
  #remove the original variable
  D_long$dep=NULL
  
  #generate the outcome (teacher ratings at waves 2,3 and 4)
  e_trating <- rnorm(dim(D_long)[1],0,sd_tratingL1)
  D_long$t_rating <- gamma0+gamma1*D_long$prev_dep+gamma2*D_long$c_age+gamma3*D_long$c_sex+gamma4*D_long$t_score.W1+gamma5*D_long$c_ses+
    gamma6*D_long$time+D_long$L3_RE_trating+ D_long$L2_RE_trating+e_trating
  
  #generate the auxiliary variable (SDQ at waves 1,2 and 3)
  e_sdq <- rnorm(dim(D_long)[1],0,sd_sdqL1)
  D_long$prev_sdq <- beta0+beta1*D_long$prev_dep+beta2*D_long$time+D_long$L3_RE_sdq+ D_long$L2_RE_sdq+e_sdq
  
  #drop variables
  D_long <- D_long %>% select(!contains("L3")& !contains("L2"))
  
  ##########MIssing data generation######################
  
  #reshape to wide format
  D_wide <- reshape(D_long,v.names=c("school","t_rating","prev_sdq","prev_dep"),timevar = "time",idvar="c_id",direction= "wide")
  
  #generate missing data in prev_dep 2,3 and 4 (probability of response, 1=observed, 0=missing) 
 
  if (waveno==4){
  D_wide$r_prevdep.2 <- as.numeric(runif(N,0,1)< inv.logit(iota01+iota1*D_wide$t_rating.2+iota2*D_wide$prev_sdq.2))
  D_wide$r_prevdep.3 <- as.numeric(runif(N,0,1)< inv.logit(iota02+iota1*D_wide$t_rating.3+iota2*D_wide$prev_sdq.3))
  D_wide$r_prevdep.4 <- as.numeric(runif(N,0,1)< inv.logit(iota03+iota1*D_wide$t_rating.4+iota2*D_wide$prev_sdq.4))
  
  
  #check missing data proportions
  1-(sum(D_wide$r_prevdep.2)/nrow(D_wide))
  1-(sum(D_wide$r_prevdep.3)/nrow(D_wide))
  1-(sum(D_wide$r_prevdep.4)/nrow(D_wide))
  } else if(waveno==8){
    
    #generate missing data in prev_dep 2-8 (probability of response, 1=observed, 0=missing) 
    D_wide$r_prevdep.2 <- as.numeric(runif(1200,0,1)< inv.logit(iota0+iota1*D_wide$t_rating.2+iota2*D_wide$prev_sdq.2))
    D_wide$r_prevdep.3 <- as.numeric(runif(1200,0,1)< inv.logit(iota0+iota1*D_wide$t_rating.3+iota2*D_wide$prev_sdq.3))
    D_wide$r_prevdep.4 <- as.numeric(runif(1200,0,1)< inv.logit(iota0+iota1*D_wide$t_rating.4+iota2*D_wide$prev_sdq.4))
    D_wide$r_prevdep.5 <- as.numeric(runif(1200,0,1)< inv.logit(iota0+iota1*D_wide$t_rating.5+iota2*D_wide$prev_sdq.5))
    D_wide$r_prevdep.6 <- as.numeric(runif(1200,0,1)< inv.logit(iota0+iota1*D_wide$t_rating.6+iota2*D_wide$prev_sdq.6))
    D_wide$r_prevdep.7 <- as.numeric(runif(1200,0,1)< inv.logit(iota0+iota1*D_wide$t_rating.7+iota2*D_wide$prev_sdq.7))
    D_wide$r_prevdep.8 <- as.numeric(runif(1200,0,1)< inv.logit(iota0+iota1*D_wide$t_rating.8+iota2*D_wide$prev_sdq.8))
    
    #check missing data proportions
    1-(sum(D_wide$r_prevdep.2)/nrow(D_wide))
    1-(sum(D_wide$r_prevdep.3)/nrow(D_wide))
    1-(sum(D_wide$r_prevdep.4)/nrow(D_wide))
    1-(sum(D_wide$r_prevdep.5)/nrow(D_wide))
    1-(sum(D_wide$r_prevdep.6)/nrow(D_wide))
    1-(sum(D_wide$r_prevdep.7)/nrow(D_wide))
    
  }
  
  D_wide <- D_wide %>% select(!c("t_rating.1","prev_sdq.1","prev_dep.1"))
  
  colnames(D_wide)[colnames(D_wide)=="school.1"] <- "school.W1"
  
  
  #reshape to long format
  
  if (waveno==4){
  D_long <- reshape(D_wide,varying =list(c("school.2","school.3","school.4"),
                                         c("prev_dep.2","prev_dep.3","prev_dep.4"),
                                         c("r_prevdep.2","r_prevdep.3","r_prevdep.4"),
                                         c("t_rating.2","t_rating.3","t_rating.4"),
                                         c("prev_sdq.2","prev_sdq.3","prev_sdq.4")),idvar="c_id", 
                    v.names=c("school","prev_dep","r_prevdep","t_rating","prev_sdq"), times=c(2,3,4),direction= "long")
  }else if(waveno==8){
    D_long <- reshape(D_wide,varying =list(c("school.2","school.3","school.4","school.5","school.6","school.7","school.8"),
                                           c("prev_dep.2","prev_dep.3","prev_dep.4","prev_dep.5","prev_dep.6","prev_dep.7","prev_dep.8"),
                                           c("r_prevdep.2","r_prevdep.3","r_prevdep.4","r_prevdep.5","r_prevdep.6","r_prevdep.7","r_prevdep.8"),
                                           c("t_rating.2","t_rating.3","t_rating.4","t_rating.5","t_rating.6","t_rating.7","t_rating.8"),
                                           c("prev_sdq.2","prev_sdq.3","prev_sdq.4","prev_sdq.5","prev_sdq.6","prev_sdq.7","prev_sdq.8")),idvar="c_id", 
                      v.names=c("school","prev_dep","r_prevdep","t_rating","prev_sdq"), times=c(2,3,4,5,6,7,8),direction= "long")
  }
  
  #assign NA values 
  D_long$prev_dep <- ifelse(D_long$r_prevdep==0,NA,D_long$prev_dep)
  
  
  D_long$r_prevdep=NULL
  #Export data
  write.csv(D_long,paste0("data",i,".csv"),row.names = F)
}

}

###########################################################################################################
#Generate the data for the different simulation scenarios
############################################################################################################

#* NOTE: Set the working directory to where data to be stored 

############### MAR-CATS- adding 10 clusters ###################################################

#MAR- CATS/ adding 10 clusters at each wave/ High-high ICC combination
simul(279109,1000,40,30,"MAR1","High","high",10,4)
  
#MAR- CATS/ adding 10 clusters at each wave/ High-low ICC combination
simul(4680,1000,40,30,"MAR1","High","low",10,4)

#MAR- CATS/ adding 10 clusters at each wave/ Low-high ICC combination
simul(627190,1000,40,30,"MAR1","Low","high",10,4)

#MAR- CATS/ adding 10 clusters at each wave/ Low-low ICC combination
simul(62517,1000,40,30,"MAR1","Low","low",10,4)

############### MAR-CATS- adding 50 clusters ###################################################

#MAR- CATS/ adding 50 clusters at each wave/ High-high ICC combination
simul(3682500,1000,40,30,"MAR1","High","high",50,4)

#MAR- CATS/ adding 50 clusters at each wave/ High-low ICC combination
simul(527819,1000,40,30,"MAR1","High","low",50,4)

#MAR- CATS/ adding 50 clusters at each wave/ Low-high ICC combination
simul(112890,1000,40,30,"MAR1","Low","high",50,4)

#MAR- CATS/ adding 50 clusters at each wave/ Low-low ICC combination
simul(11902,1000,40,30,"MAR1","Low","low",50,4)

############### MAR-inflated- adding 10 clusters ###################################################

#MAR-inflated/ adding 10 clusters at each wave/ High-high ICC combination
simul(226109,1000,40,30,"MAR2","High","high",10,4)

#MAR-inflated/ adding 10 clusters at each wave/ High-low ICC combination
simul(726199,1000,40,30,"MAR2","High","low",10,4)

#MAR-inflated/ adding 10 clusters at each wave/ Low-high ICC combination
simul(2518018,1000,40,30,"MAR2","Low","high",10,4)

#MAR-inflated/ adding 10 clusters at each wave/ Low-low ICC combination
simul(179291,1000,40,30,"MAR2","Low","low",10,4)

#MAR-inflated/ adding 50 clusters at each wave/ High-high ICC combination
simul(7341680,1000,40,30,"MAR2","High","high",50,4)

#MAR-inflated/ adding 50 clusters at each wave/ High-low ICC combination
simul(425719,1000,40,30,"MAR2","High","low",50,4)

#MAR-inflated/ adding 50 clusters at each wave/ Low-high ICC combination
simul(268103,1000,40,30,"MAR2","Low","high",50,4)

#MAR-inflated/ adding 50 clusters at each wave/ Low-low ICC combination
simul(80298,1000,40,30,"MAR2","Low","low",50,4)

############### MAR-inflated- small sample ###################################################

#MAR-inflated/ small sample/ High-high ICC combination
simul(3781010,1000,40,7.5,"MAR2","High","high",10,4)

#MAR-inflated/ small sample/ High-low ICC combination
simul(92010837,1000,40,7.5,"MAR2","High","low",10,4)

#MAR-inflated/ small sample/ Low-high ICC combination
simul(6728903,1000,40,7.5,"MAR2","Low","high",10,4)

#MAR-inflated/ small sample/ Low-low ICC combination
simul(83927640,1000,40,7.5,"MAR2","Low","low",10,4)

############### MAR-inflated- higher waves ###################################################

#MAR-inflated/ higher waves/ High-high ICC combination
simul(3729222,1000,40,30,"MAR2","High","high",10,8)

#MAR-inflated/ higher waves/ High-low ICC combination
simul(6735241,1000,40,30,"MAR2","High","low",10,8)

#MAR-inflated/ higher waves/ Low-high ICC combination
simul(7382111,1000,40,30,"MAR2","Low","high",10,8)

#MAR-inflated/ higher waves/ Low-low ICC combination
simul(94021843,1000,40,30,"MAR2","Low","low",10,8)


############### MAR-inflated- constant number of higher-level clusters ##########################

#MAR-inflated/ constant number of higher-level clusters/ High-high ICC combination
simul(9674265,1,10,120,"MAR2","High","high",0,4)

#MAR-inflated/ constant number of higher-level clusters/ High-low ICC combination
simul(28101444,1000,10,120,"MAR2","High","low",0,4)

#MAR-inflated/ constant number of higher-level clusters/ Low-high ICC combination
simul(38244427,1000,10,120,"MAR2","Low","high",0,4)

#MAR-inflated/ constant number of higher-level clusters/ Low-low ICC combination
simul(474229022,1000,10,120,"MAR2","Low","low",0,4)


