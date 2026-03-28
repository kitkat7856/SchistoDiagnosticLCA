## LOAD DATA 
#### Kato Katz Data ####

loadKK <- function(dataframe, ids){
  dt <- kkdata # read data file
  
  # name different weeks
  dt$tmpN <- NA
  dt$tmpN[which(dt$date=="25/09/2017" | dt$date=="26/09/2017" | dt$date=="27/09/2017"
                | dt$date=="28/09/2017" | dt$date=="29/09/2017" | dt$date=="02/10/2017")] <- "Pre-T"
  
  dt$tmpN[which(dt$date=="23/10/2017" | dt$date=="24/10/2017" | dt$date=="25/10/2017"
                | dt$date=="26/10/2017" | dt$date=="27/10/2017")] <- "3 weeks"
  
  dt$tmpN[which(dt$date=="04/12/2017" | dt$date=="05/12/2017")] <- "9 weeks"
  
  dt$tmpN[which(dt$date=="01/03/2018" | dt$date=="05/03/2018" | dt$date=="06/03/2018"
                | dt$date=="07/03/2018" | dt$date=="08/03/2018" | dt$date=="09/03/2018")] <- "6 months"
  
  tmps <- c(na.omit(unique(dt$tmpN)))
  
  # set as array dimension (number child_ID, number timesteps, number max of repeats)
  
  KK <- array(NA,dim = c(length(ids),length(tmps),6))
  
  for (i in 1:length(ids)){ #for each child
    dts<- subset(dt, child_id==ids[i]) #temporary data frame for child ID
    KKchild<- matrix(NA, length(tmps), 6) #matrix for 6 KK measures, for each child
    for(t in 1:length(tmps)){
      kkrepeats <-c(dts$Sm_A[which(dts$tmpN==tmps[t])], dts$Sm_B[which(dts$tmpN==tmps[t])]) #KK measurements for a timepoint
      length(kkrepeats)<-6 #max number repeats is 6
      KKchild[t,]<- kkrepeats
    }
    KK[i,,] <- KKchild
  }
  return(KK)
}

#### CCA Data ####

loadCCA <- function(dataframe,ids){
  dt <- dataframe
  
  dt$cca<- NA
  dt$cca[which(dt$Poppy == "I*" | dt$Poppy=="-")]<-NA
  
  # restructure for categorical
  dt$cca[which(dt$Poppy=="T")] <- 1 #trace
  dt$cca[which(dt$Poppy==3)]<-4 #+++
  dt$cca[which(dt$Poppy==2)]<-3 #++
  dt$cca[which(dt$Poppy==1)]<-2 #+
  dt$cca[which(dt$Poppy==0.5)]<-1 #trace
  dt$cca[which(dt$Poppy==0)]<-0 #negative
  
  # name different timepoints
  dt$tmpN <- NA
  dt$tmpN[which(dt$date_on_tube=="25/09/2017" | dt$date_on_tube=="26/09/2017" | dt$date_on_tube=="27/09/2017"
                | dt$date_on_tube=="28/09/2017" | dt$date_on_tube=="29/09/2017" | dt$date_on_tube=="02/10/2017")] <- "Pre-T"
  
  dt$tmpN[which(dt$date_on_tube=="23/10/2017" | dt$date_on_tube=="24/10/2017" | dt$date_on_tube=="25/10/2017"
                | dt$date_on_tube=="26/10/2017" | dt$date_on_tube=="27/10/2017")] <- "3 weeks"
  
  dt$tmpN[which(dt$date_on_tube=="04/12/2017" | dt$date_on_tube=="05/12/2017" | dt$date_on_tube=="4/12/17" | dt$date_on_tube=="5/12/17")] <- "9 weeks"
  
  dt$tmpN[which(dt$date_on_tube=="01/03/2018" | dt$date_on_tube=="05/03/2018" | dt$date_on_tube=="06/03/2018"
                | dt$date_on_tube=="07/03/2018" | dt$date_on_tube=="08/03/2018" | dt$date_on_tube=="09/03/2018")] <- "6 months"
  dt$tmpN[which(dt$date_on_tube=="1/3/18" | dt$date_on_tube=="5/3/18" | dt$date_on_tube=="6/3/18"
                | dt$date_on_tube=="7/3/18" | dt$date_on_tube=="8/3/18" | dt$date_on_tube=="9/3/18")] <- "6 months"
  
  tmps <- c(na.omit(unique(dt$tmpN)))
  
  # set as array dimension (number child_ID, number timesteps)
  CCA <- matrix(NA,length(ids),length(tmps))
  
  for (i in 1:length(ids)){ #for each child
    dts<- subset(dt, CID==ids[i]) #temporary dataframe for child i
    ccaChild<- c() #vector of cca measurements of child i
    
    for (t in 1:length(tmps)){
      ccaChild[t]<- ifelse(is.null(dts$cca[which(dts$tmpN == tmps[t])]), NA,dts$cca[which(dts$tmpN==tmps[t])])
    }
    
    CCA[i,] <- ccaChild
  }
  return(CCA)
} 

#### CCA GScore data ####

loadGscore <- function(dataframe,ids){
  #restructure so that it is the right scoring so that a
  # score of 1 which is negative, is given a score of 0.
  dt <- dataframe
  
  dt$gscore <- dt$gscore-1
  
  # name different weeks
  dt$tmpN <- NA
  dt$tmpN[which(dt$date_on_tube=="25/09/2017" | dt$date_on_tube=="26/09/2017" | dt$date_on_tube=="27/09/2017"
                | dt$date_on_tube=="28/09/2017" | dt$date_on_tube=="29/09/2017" | dt$date_on_tube=="02/10/2017")] <- "Pre-T"
  
  dt$tmpN[which(dt$date_on_tube=="23/10/2017" | dt$date_on_tube=="24/10/2017" | dt$date_on_tube=="25/10/2017"
                | dt$date_on_tube=="26/10/2017" | dt$date_on_tube=="27/10/2017")] <- "3 weeks"
  
  dt$tmpN[which(dt$date_on_tube=="04/12/2017" | dt$date_on_tube=="05/12/2017" | dt$date_on_tube=="4/12/17" | dt$date_on_tube=="5/12/17")] <- "9 weeks"
  
  dt$tmpN[which(dt$date_on_tube=="01/03/2018" | dt$date_on_tube=="05/03/2018" | dt$date_on_tube=="06/03/2018"
                | dt$date_on_tube=="07/03/2018" | dt$date_on_tube=="08/03/2018" | dt$date_on_tube=="09/03/2018")] <- "6 months"
  dt$tmpN[which(dt$date_on_tube=="1/3/18" | dt$date_on_tube=="5/3/18" | dt$date_on_tube=="6/3/18"
                | dt$date_on_tube=="7/3/18" | dt$date_on_tube=="8/3/18" | dt$date_on_tube=="9/3/18")] <- "6 months"
  
  tmps <- c(na.omit(unique(dt$tmpN)))
  
  # set as array dimension (number child_ID, number timesteps)
  gscore <- array(NA,dim = c(length(ids),length(tmps)))
  
  for (i in 1:length(ids)){ #for each child
    dts<- subset(dt, CID==ids[i]) #temporary dataframe for child i
    gscoreChild<- c() #vector of gscore measurements of child i
    
    for (t in 1:length(tmps)){ #for each timepoint
      gscoreChild[t]<- ifelse(is.null(dts$gscore[which(dts$tmpN == tmps[t])]), NA, as.numeric(dts$gscore[which(dts$tmpN==tmps[t])]))
    }
    
    gscore[i,] <- gscoreChild
  }
  
  return(gscore)
} 

####  qpcr ####
loadCT<- function(dataframe,ids){
  #stool qpcr CT data into matrix
  require(dplyr)
  dt<- dataframe
  dt$CT[which(dt$CT == "Undetermined")]<- 45
  dt$CT<- as.numeric(dt$CT)
  
  dts<- dt %>%
    group_by(Sample.Name) %>%
    summarise(meanCT = mean(CT, na.rm=T)) %>%
    mutate(tmp=sub("_.*", "",Sample.Name), CID=sub(".*_", "",Sample.Name))
  
  #matrix for CT, dimensions ID, timepoints
  timepoints<- c("preT", "3wk", "9wk", "6mo")
  
  CT<- matrix(NA, length(ids), length(timepoints))  
  
  for (i in 1:length(ids)){
    CTchild<- subset(dts, dts$CID== ids[i])
    CTtmp<- c()
    
    for(t in 1:length(timepoints)){
      CTtmp[t]<- ifelse(is.null(CTchild$meanCT[which(CTchild$tmp == timepoints[t])]), NA, CTchild$meanCT[which(CTchild$tmp == timepoints[t])])
    }
    
    CT[i,]<- CTtmp
  }
  
  return(CT)
}

loadDBS<- function(dataframe, ids, col){
  #stool qpcr CT data into matrix
  require(dplyr)
  dt<- dataframe
  dt$CT[which(dt$CT == "Undetermined")]<- 45
  dt$CT<- as.numeric(dt$CT)
  
  dts<- dt %>%
    group_by(Sample.Name) %>%
    summarise(meanCT=mean(CT, na.rm=T)) %>%
    mutate(tmp= sub("_.*", "",Sample.Name),CID=sub(".*_", "",Sample.Name))
  
  dts$binary<- NA
  dts$binary[dts$meanCT>=38]<-0
  dts$binary[dts$meanCT<38]<-1
  
  #matrix for CT, dimensions ID, timepoints
  timepoints<- c("preT", "3wk", "9wk", "6mo")
  
  DBS<- matrix(NA, length(ids), length(timepoints))  
  
  for (i in 1:length(ids)){
    DBSchild<- subset(dts, dts$CID== ids[i])
    DBStmp<- matrix(NA,4,2) #1 mean ct, #2 mean binary
    
    for(t in 1:length(timepoints)){
      DBStmp[t,1]<- ifelse(is.null(DBSchild$meanCT[which(DBSchild$tmp == timepoints[t])]), NA, DBSchild$meanCT[which(DBSchild$tmp == timepoints[t])])
      DBStmp[t,2]<- ifelse(is.null(DBSchild$binary[which(DBSchild$tmp == timepoints[t])]), NA, DBSchild$binary[which(DBSchild$tmp == timepoints[t])])
    }
    
    DBS[i,]<- DBStmp[,col]
  }
  
  return(DBS)
}

loadHatch<- function(dataframe, ids,col){
  require(stringr)
  require(dplyr)
  
  #set data file
  dt<- dataframe
  dt$CID<- as.numeric(str_sub(dt$id_full, 6, 9))
  
  #set tmps
  dt$tmpN <- NA
  dt$tmpN[which(dt$date_collected=="25/09/2017" | dt$date_collected=="26/09/2017" | dt$date_collected=="27/09/2017"
                | dt$date_collected=="28/09/2017" | dt$date_collected=="29/09/2017" | dt$date_collected=="02/10/2017")] <- "Pre-T"
  
  dt$tmpN[which(dt$date_collected=="23/10/2017" | dt$date_collected=="24/10/2017" | dt$date_collected=="25/10/2017"
                | dt$date_collected=="26/10/2017" | dt$date_collected=="27/10/2017")] <- "3 weeks"
  
  dt$tmpN[which(dt$date_collected=="04/12/2017" | dt$date_collected=="05/12/2017")] <- "9 weeks"
  
  dt$tmpN[which(dt$date_collected=="05/03/2018" | dt$date_collected=="06/03/2018"| dt$date_collected=="07/03/2018"
                | dt$date_collected=="08/03/2018" | dt$date_collected=="09/03/2018")] <- "6 months"
  
  tmps <- c(na.omit(unique(dt$tmpN)))
  
  
  #make matrix
  HATCH<- array(NA, dim=c(length(ids), length(tmps),2)) 
  
  for (i in 1:length(ids)){
    HATCHchild<- subset(dt, dt$CID== ids[i])
    HATCHtmp<- c()
    for(t in 1:length(tmps)){
      if(tmps[t] %in% HATCHchild$tmpN){
        HATCHtmp[t]<- sum(as.numeric(HATCHchild$mira_stored[which(HATCHchild$tmpN == tmps[t])]), na.rm=T)
      } else{
        HATCHtmp[t]<-NA
      }
    }
    HATCH[i,,1]<-HATCHtmp
  }
  HATCH[,,2][HATCH[,,1]==0]<-0
  HATCH[,,2][HATCH[,,1]>0]<-1
  HATCH[,,2][is.na(HATCH[,,1])]<-NA
  
  return(as.matrix(HATCH[,,col]))
}

#### MODEL ######

#probDBS[n,t,1] <- pneg
#probDBS[n,t,2] <- 1/(1+exp(-k2*(lambda[n,t]-int2)))
#probDBS[n,t,3] <- 1/(1+exp(-k2*(lambda[n,t]-int2)))
#DBS[n,t] ~ dbern(probDBS[n,t,status[n,t]])

### model with random walk connecting timepoints, gscore linked to lambda
eggworm <- function (N, Ti, KK, Gscore, CT, Hatch) {
  ## Set seed ##
  ##this is so that they start on different values ##
  .RNG.seed <- function(chain)
    return( switch(chain, "1"= 1, "2"= 2) )
  .RNG.name <- function(chain)
    return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )
  sigma<- c(rep(NA,3),0.01)
  
  ## Initialize Status ##
  status <- matrix(3,nrow=N, ncol=Ti)
  
  m <- "model {

   # Prior prevalence # # 1- uninfected, 2- infected, not sheddding, 3- infected and shedding
   # proportion of population at each tmp of each infection catagory 
   
    for(t in 1:Ti){
      prob[t,1:3] ~ ddirch(c(1,1,1)) #dirichlet distribution - like multiple beta 
    }
     
    # Prior KKs 
    rtnb ~ dgamma(286.3235, 223.4650)
  
    #lambda (egg infection intensity) sh/rt
    sh ~ dgamma(83.88592,125.70331)
    rt1 ~ dbeta(47.13542,633.08366)
    rt <- rt1/(1-rt1)
    
    #stool qpcr 
    a_ct ~ dgamma(0.001,0.001)
    b_ct ~ dnorm(0,0.001)
    
    #gscore cca  
    k~ dgamma(0.001,0.001)
    int ~ dnorm(0,0.001)
    negcca~dunif(0,1)
    
    #hatch binary
    #k2~ dnorm(0,0.001)
    #int2 ~ dnorm(0,0.001)
    phatch~dbeta(10,2)
    
    #random walk
    sigma[1] <- 0 
    sigma[2] <- sigma[4]*3/15
    sigma[3] <- sigma[4]*6/15
    sigma[4] ~ dgamma(0.001, 0.001) 
    
    
    #MODEL
    for(n in 1:N){
    
    ###baseline values T=1
      tEgg[n,1,1]<-0
      tEgg[n,1,2]<-0
      tEgg[n,1,3] ~ dgamma(sh,rt)
    
      for (t in 2:Ti){
        tEgg[n,t,1]<-0
        tEgg[n,t,2]<-0
        tEgg[n,t,3] ~ dnorm(tEgg[n,t-1,status[n,t-1]], sigma[t])T(0,)
      }
  
      for (t in 1:Ti){
        status[n,t] ~ dcat(prob[t,1:3]) #categorical distribution - probability of status 1 is prob[1] , status 2 is prob[2], status 3 is prob[3]
        
        lambda[n,t]<- tEgg[n,t,status[n,t]]
        for(r in 1:6){
          KK[n,t,r] ~ dnegbin(rtnb/(lambda[n,t]+rtnb),rtnb)
        }
        
        #relate ct to egg lambda- so only those with status 3(shedding) are positive, 1&2 are negative for egg
        
        meanct[n,t]<- -a_ct*log(lambda[n,t]+0.01)+b_ct ##so if egg = 0 ct is intercept (with some error)
        CT[n,t] ~ dnorm(meanct[n,t],0.08)
        
        probCCA[n,t,1] <- negcca
        probCCA[n,t,2] <- 9/(1+exp(-k*(lambda[n,t]-int)))
        probCCA[n,t,3] <- 9/(1+exp(-k*(lambda[n,t]-int)))
        
        Gscore[n,t] ~ dnorm(probCCA[n,t,status[n,t]], 1.093606)
         
        probHatch[n,t,1]<- 0
        probHatch[n,t,2]<- 0
        probHatch[n,t,3]<- phatch
        
        #probHatch[n,t,1]<- 0
        #probHatch[n,t,2]<- 0
        #probHatch[n,t,3]<- 1-exp(-exp(k2*(lambda[n,t])+ int2))
        
        Hatch[n,t] ~ dbern(probHatch[n,t,status[n,t]])
      } 
    }
    
 
  #inits# .RNG.seed, .RNG.name,status
  #data# N, Ti, KK, Gscore, CT, Hatch
  #monitor# prob, sigma, rtnb, sh, rt, k, int, negcca, a_ct, b_ct,phatch, status
  
}"
  
  # Run model #
  Results <- run.jags(m, burnin=10000, sample=20000, n.chains=2, jags.refresh = 1, method = 'parallel',
                      plots = F, silent.jags = F, summarise = T)
  return(Results)
} 

###### ROC FUNCTIONS ######

createROC<- function(diagnostic, timepoint, binary_status){
  if(length(unique(as.vector(as.matrix(binary_status))))!=2){
    stop("Status must be binary")
  }
  
  a<- 1 + 210*(timepoint-1)
  b<- 210*timepoint
  
  l<- binary_status[,a:b]
  
  prediction<- rep(list(na.omit(diagnostic[,timepoint])),500)
  label<- l[,which(!is.na(diagnostic[,timepoint]))]
  label_list<- as.list(as.data.frame(t(label)))
  
  pred<-prediction(prediction,label_list)
  roc<- performance(pred,"tpr","fpr")
  
  return(roc)
}
createSS<- function(diagnostic, timepoint, binary_status){
  if(length(unique(as.vector(as.matrix(binary_status))))!=2){
    stop("Status must be binary")
  }
  
  a<- 1 + 210*(timepoint-1)
  b<- 210*timepoint
  
  l<- binary_status[,a:b]
  
  prediction<- rep(list(na.omit(diagnostic[,timepoint])),500)
  label<- l[,which(!is.na(diagnostic[,timepoint]))]
  label_list<- as.list(as.data.frame(t(label)))
  
  pred<-prediction(prediction,label_list)
  ss<- performance(pred, "sens", "spec")
  
  return(ss)
}
getoptcut<- function(ss, samples){
  optse<- c()
  optsp<- c()
  optcut2<-c()
  for(i in 1:samples){
    cutpoints<- c(unlist(ss@alpha.values[[i]]))
    sensitivity<- c(unlist(ss@y.values[[i]]))
    specificity<- c(unlist(ss@x.values[[i]]))
    
    d<- sqrt((1-sensitivity)^2 + (1-specificity)^2)
    optcut2[i]<- cutpoints[which.min(d)]
    optse[i]<- sensitivity[which(cutpoints == optcut2[i])]
    optsp[i]<- specificity[which(cutpoints == optcut2[i])]
    
  }
  d<- data.frame(cut2= mean(optcut2), se= mean(optse), sp=mean(optsp), sdse=sd(optse), sdsp=sd(optsp))
  print(d)
}


